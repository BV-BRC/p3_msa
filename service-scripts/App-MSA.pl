#The MSA application with variance analysis.

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use P3DataAPI;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use LWP::UserAgent;
use JSON::XS;
use JSON;
use IPC::Run qw(run);
use Cwd;
use Clone;
use URI::Escape;

my $script = Bio::KBase::AppService::AppScript->new(\&process_fasta);
my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;

my %aacode = (
TTT => "F", TTC => "F", TTA => "L", TTG => "L",
TCT => "S", TCC => "S", TCA => "S", TCG => "S",
TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
TGT => "C", TGC => "C", TGA => "*", TGG => "W",
CTT => "L", CTC => "L", CTA => "L", CTG => "L",
CCT => "P", CCC => "P", CCA => "P", CCG => "P",
CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
CGT => "R", CGC => "R", CGA => "R", CGG => "R",
ATT => "I", ATC => "I", ATA => "I", ATG => "M",
ACT => "T", ACC => "T", ACA => "T", ACG => "T",
AAT => "N", AAC => "N", AAA => "K", AAG => "K",
AGT => "S", AGC => "S", AGA => "R", AGG => "R",
GTT => "V", GTC => "V", GTA => "V", GTG => "V",
GCT => "A", GCC => "A", GCA => "A", GCG => "A",
GAT => "D", GAC => "D", GAA => "E", GAG => "E",
GGT => "G", GGC => "G", GGA => "G", GGG => "G",
);

my $THRESHOLD = 0.75;

my $rc = $script->run(\@ARGV);

exit $rc;

sub process_fasta
{
    my($app, $app_def, $raw_params, $params) = @_;
    print "Proc MSA Var ", Dumper($app_def, $raw_params, $params);
    my $token = $app->token();
    my $data_api_module = P3DataAPI->new($data_api, $token);
    my $output_folder = $app->result_folder();
    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    #
    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";
    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";
    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);
    #
    # Read parameters and discover input files that need to be staged.
    #
    # Make a clone so we can maintain a list of refs to the paths to be
    # rewritten.
    #
    my %in_files;
    my $params_to_app = Clone::clone($params);
    #
    # Count the number of files.
    #
    my $file_count = 0;
    if (exists($params_to_app->{fasta_files})) {
        $file_count = $file_count + scalar(@{$params_to_app->{fasta_files}});
    }
    if (exists($params_to_app->{feature_groups})) {
        $file_count = $file_count + scalar(@{$params_to_app->{feature_groups}});
    }
    if (length($params_to_app->{fasta_keyboard_input}) >= 1) {
        $file_count = $file_count + 1;
    }
    say STDERR "Number of files: $file_count.";
    my $prefix = $params_to_app->{output_file};
    #
    # Determine if the data is represented as DNA or protein.
    #
    my $dna = 1; # Use the DNA alphabet.
    my $in_type = "feature_dna_fasta";
    if (substr($params_to_app->{alphabet}, 0, 1) eq "p") {
    	$dna = 0; # Use the amino acid, protein alphabet.
	    $in_type = "feature_protein_fasta";
    }
    my @genome_groups;
    if (exists($params_to_app->{select_genomegroup})) {
        $dna = 1;
        $in_type = "feature_dna_fasta";
        $file_count = $file_count + scalar(@{$params_to_app->{select_genomegroup}});
        for my $group_path (@{$params_to_app->{select_genomegroup}}) {
            say STDERR "Getting genome group: $group_path";
            my $genome_id_exclude = "";
            if ($params_to_app->{ref_type} eq "genome_id") {
                $genome_id_exclude = $params_to_app->{ref_string};
                $genome_id_exclude =~ s/^\s+|\s+$//g;
            }
            push @genome_groups, get_genome_group_file($data_api_module,
                                                       $group_path,
                                                       $stage_dir,
                                                       $genome_id_exclude
                                                       );
        }
    }
    #
    # Write files to the staging directory.
    #
    my @to_stage;
    my $aligned_exists = 0;
    my $mixed = 0;
    for my $read_tuple (@{$params_to_app->{fasta_files}}) {
        for my $read_name (keys %{$read_tuple}) {
            if($read_name eq "file") {
                my $nameref = \$read_tuple->{$read_name};
                $in_files{$$nameref} = $nameref;
                push(@to_stage, $$nameref);
            }
            else {
                if(index($read_tuple->{$read_name}, "aligned") != -1) {
                    $aligned_exists = 1;
                }
                if (index($read_tuple->{$read_name}, "protein") != -1) {
                    $mixed = 1;
                }
            }
        }
    }
    if ($mixed == 1) {
        $dna = 0;
        $in_type = "feature_protein_fasta";
    }
    my $staged = {};
    if (@to_stage)
    {
        warn Dumper(\%in_files, \@to_stage);
        $staged = $app->stage_in(\@to_stage, $stage_dir, 1);
        while (my($orig, $staged_file) = each %$staged)
        {
            my $path_ref = $in_files{$orig};
            $$path_ref = $staged_file;
        }
    }
    #
    # Download feature groups in a file.
    #
    my $ofile = "$stage_dir/feature_groups.fasta";
    open(F, ">$ofile") or die "Could not open $ofile";
    my $feature_id_exclude = "";
    if ($params_to_app->{ref_type} eq "feature_id") {
        $feature_id_exclude = $params_to_app->{ref_string};
        $feature_id_exclude =~ s/^\s+|\s+$//g;
    }
    for my $feature_name (@{$params_to_app->{feature_groups}}) {
	    my $ids = $data_api_module->retrieve_patricids_from_feature_group($feature_name);
        my @ids_new = ();
        for my $id (@$ids){
	        if ($id ne $feature_id_exclude) {
                push(@ids_new, $id);
            }
        }
	    my $seq = "";
	    if ($dna) {
		    $seq = $data_api_module->retrieve_nucleotide_feature_sequence(\@ids_new);
	    } else {
		    $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids_new);
	    }
	    for my $id (@ids_new) {
		    my $out = ">$id\n" . $seq->{$id} . "\n";
    		    print F $out;
	    }
    }
    if (exists($params_to_app->{feature_groups})) {
    	push @{ $params_to_app->{fasta_files} }, {"file" => $ofile, "type" => $in_type};
    }
    close(F);
    #
    # Put keyboard input into a file.
    #
    my $text_input_file = "$stage_dir/fasta_keyboard_input.fasta";
    if ((not $dna) && $params_to_app->{fasta_keyboard_input} && (not is_aa($params_to_app->{fasta_keyboard_input}, 0))) {
        convert_aa_file($params_to_app->{fasta_keyboard_input}, $text_input_file, 0);
    } else {
        open(FH, '>', $text_input_file) or die "Cannot open $text_input_file: $!";
        print FH $params_to_app->{fasta_keyboard_input};
        close(FH);
    }
    push @{ $params_to_app->{fasta_files} }, {"file" => $text_input_file, "type" => $in_type};
    #
    # Add genome groups.
    #
    for my $group_file (@genome_groups) {
        push @{ $params_to_app->{fasta_files} }, {"file" => $group_file,
                                                  "type" => $in_type};
    }
    #
    # Check if there is a protein file even if we are supposed to be doing DNA. Set to protein if so.
    # This ensures that the afa file type is set correctly for the MSA viewer.
    #
    if ($dna && (scalar(@genome_groups) <= 0)) {
        for my $read_tuple (@{$params_to_app->{fasta_files}}) {
            if (is_aa($read_tuple->{file})) {
                say STDERR "Expecting DNA files, but the file is aa: " . $read_tuple->{file};
                $dna = 0;
                last;
            }
        }
    }
    say STDERR "Alignment already present?: $aligned_exists";
    say STDERR "Using DNA?: $dna"; #  Protein files exist: $mixed";
    say STDERR "Mixed?: $mixed";
    #
    # Combine all files into one input.fasta file.
    # Put the reference sequence first if present.
    #
    my $work_fasta = "$work_dir/input.fasta";
    open(IN, '>', $work_fasta) or die "Cannot open $work_fasta: $!";
    my $ref_string = $params_to_app->{ref_string};
    $ref_string =~ s/^\s+|\s+$//g;
    if ($params_to_app->{ref_type} eq "string") {
        print IN $ref_string . "\n";
    } elsif ($params_to_app->{ref_type} eq "feature_id") {
        my @ids = ($ref_string);
        my $seq = "";
        if ($dna == 1) {
            $seq = $data_api_module->retrieve_nucleotide_feature_sequence(\@ids);
        } else {
            $seq = $data_api_module->retrieve_protein_feature_sequence(\@ids);
        }
        print IN ">" . $ref_string . "\n" . $seq->{$ref_string} . "\n" or die "Could not write reference sequence.";
    } elsif ($params_to_app->{ref_type} eq "genome_id") {
        my @ids = ($ref_string);
        my $reference_genome_file = get_genome_seqs($data_api_module, \@ids, \*IN, $stage_dir, "");
    }
    for my $read_tuple (@{$params_to_app->{fasta_files}}) {
    	my $filename = $read_tuple->{file};
        my $convert = 0;
        if ((index($read_tuple->{type}, "dna") != -1) && (not $dna) && not is_aa($filename)) {
            $convert = 1;
        }
        open my $fh, '<', $filename or die "Cannot open $filename: $!";
        my $seq_line = "";
        while ( my $line = <$fh> ) {
            # chomp($line); # remove newlines
            # $line =~ s/#.*//; # remove comments
            # $line =~ s/;.*//; # remove comments
            # $line =~ s/^\s+//;  # remove leading whitespace
            # $line =~ s/\s+$//; # remove trailing whitespace
            my $print_me = 1;
            next if(length($line) <= 0);
            if ($aligned_exists && $file_count > 1 && substr($line, 0, 1) ne ">") {
                $line =~ tr/-_.~*//d; # Remove indels from alignments if other files are present.
            }
            if ($convert && substr($line, 0, 1) ne ">") {
                chomp($line);
                $seq_line = $seq_line . $line;
                $print_me = 0;
            } elsif ($convert) {
                if ($seq_line) {
                    print IN convert_aa_line($seq_line) . "\n";
                }
                $seq_line = "";
            }
            if ($print_me) {
                print IN $line;
            }
        }
        if ($seq_line) {
            print IN convert_aa_line($seq_line) . "\n";
        }
        close($fh);
    }
    close(IN);
    #
    # Run the multiple sequence aligner.
    #
    my $recipe = lc($params_to_app->{aligner});
    if ($aligned_exists && $file_count == 1) {
        rename "$work_dir/input.fasta", "$work_dir/output.afa";
    }
    elsif ($recipe eq "muscle") {
        # An aln file only displayed 4 out of 8 sequences in MView and JalView for some reason. I could not find anything wrong with the format. Removing clustal w format.
    	print STDOUT "Running MUSCLE.\n";
        my @muscle_version = ("muscle", "-version");
        my @muscle_cmd =  ("muscle", "-quiet", "-in", "$work_dir/input.fasta", "-fastaout", "$work_dir/output.afa"); # , "-clwstrict", "-clwout", "$work_dir/$prefix.aln");
        my $string_cmd = join(" ", @muscle_cmd);
        run_cmd(\@muscle_version);
        run(\@muscle_version, "1>>", "$work_dir/muscle.job.log");
        open(MUSCLE_LOG, '>>', "$work_dir/muscle.job.log") or die $!;
        print MUSCLE_LOG "$string_cmd\n";
        close(MUSCLE_LOG) or die $!;
        print STDOUT "$string_cmd\n";
        run_cmd(\@muscle_cmd);
        print STDOUT "Finished MUSCLE.\n";
    }
    elsif ($recipe eq "progressivemauve" or $recipe eq "mauve") {
        # The file format is different.
        my $mauve_out = "$work_dir/$prefix.xmfa";
        my @mauve_cmd = ("progressiveMauve", "--output=$mauve_out", "$work_dir/input.fasta");
        my $string_cmd = join(" ", @mauve_cmd);
        print STDOUT "Running progressiveMauve.\n";
        print STDOUT "$string_cmd\n";
        my $ok = run(\@mauve_cmd, "&>", "$work_dir/$prefix.mauve.log");
        if (!$ok) {
            die "progressiveMauve command failed.\n";
        }
        print STDOUT "Finished progressiveMauve.\n";
        my @id_arr = ();
        open(INF, '<', "$work_dir/input.fasta") or die "Couldn't open file $work_dir/input.fasta. $!";
        while (my $line = <INF>) {
            if (substr($line, 0, 1) eq ">") {
                push(@id_arr, $line);
            }
        }
        close(INF) or die $!;
        open(INF, '<', $mauve_out) or die "Couldn't open file $mauve_out. $!";
        open(OUTF, '>', "$work_dir/output.afa") or die "Couldn't open file $work_dir/output.afa. $!";
        my $printme = 0;
        while (my $line = <INF>) {
            if (substr($line, 0, 1) eq ">") {
                print OUTF $id_arr[$printme];
                $printme = $printme + 1;
            } elsif (substr($line, 0, 1) eq "=") {
                last;
            } elsif ($printme) {
                print OUTF $line;
            }
        }
        close(INF) or die $!;
        close(OUTF) or die $!;
    }
    elsif ($recipe eq "mafft") {
        print STDOUT "Running mafft.\n";
        my @mafft_cmd = ("mafft", "--auto", "--preservecase", "$work_dir/input.fasta");
        my $string_cmd = join(" ", @mafft_cmd);
        my @mafft_version = ("mafft", "--version");
        run(\@mafft_version, "2>>", "$work_dir/mafft.job.log");
        run_cmd(\@mafft_version);
        open(MAFFT_LOG, '>>', "$work_dir/mafft.job.log") or die $!;
        print MAFFT_LOG "$string_cmd\n";
        close(MAFFT_LOG) or die $!;
        print STDERR "$string_cmd\n";
        my $ok = run(\@mafft_cmd, "1>", "$work_dir/output.afa", "2>>", "$work_dir/mafft.job.log");
        if (!$ok) {
            die "Mafft command failed.\n";
        }
        print STDOUT "Finished mafft.\n"
    }
    else {
        die "Recipe not found: $recipe\n";
    }
    rename "$work_dir/output.afa", "$work_dir/$prefix.afa";
    open(my $prefix_afa, "<",  "$work_dir/$prefix.afa") or die "Couldn't open $work_dir/$prefix.afa";
    open(my $output_afa, ">", "$work_dir/output.afa") or die "Couldn't open $work_dir/output.afa";
    while ( my $line = <$prefix_afa> ) {
        next if(length($line) <= 0);
        if (substr($line, 0, 1) eq ">") {
            print $output_afa $line;
        } else {
            print $output_afa uc($line);
        }
    }
    close($prefix_afa);
    close($output_afa);
    #
    # Run the SNP analysis.
    # Requires the alignment file to be named 'output.afa'
    #
    my @cmd = ("snp_analysis", "-r", "$work_dir", "-x");
    if ($dna) {
    	push @cmd, "-n";
    }
    run_cmd(\@cmd);
    unlink("$work_dir/output.afa") or warn "Unable to unlink $output_afa: $!";
    print STDOUT "Completed SNP analysis.\n";
    # my $file_str = "$work_dir/cons.fasta";
    # my $cons_out = "$work_dir/$prefix.consensus.fasta";
    # open(INC, '<', "$work_dir/cons.fasta") or die "Couldn't open file $work_dir/cons.fasta. $!";
    # open(CONS, '>', "$work_dir/$prefix.consensus.fasta") or die "Couldn't open file $work_dir/$prefix.consensus.fasta. $!";
    # my $line_num = 0;
    # while (my $line = <INC>) {
    #     if ($line_num == 0) {
    #         my @spl = split('|', $line);
    #         my $head_str = $spl[0]."|"."num_sequences:".$spl[1]."|"."alignment_length:".$spl[2]."|"."consensus_length:".$spl[3];
    #         print CONS $head_str;
    #     } else {
    #         print CONS $line;
    #     }
    #     print STDERR $line;
    #     $line_num = $line_num + 1;
    # }
    # close INC or die $!;
    # close CONS or die $!;
    # unlink(\$file_str) or warn "Unable to unlink $file_str: $!";
    rename "$work_dir/cons.fasta", "$work_dir/$prefix.consensus.fasta";

    my @cmd = ("convert_seq_files", "$work_dir/$prefix.afa", "$work_dir/$prefix");
    if ($dna) {
    	push @cmd, "DNA";
    } else {
        push @cmd, "protein";
    }
    run_cmd(\@cmd);
    #
    # Create figures.
    #
    @cmd = ("snp_analysis_figure", "$work_dir/foma.table", "$work_dir/$prefix");
    run_cmd(\@cmd);
    print STDOUT "Completed figure creation.\n";
    #
    # Add a position relative to the reference sequence if it exists.
    #
    if ($params_to_app->{ref_type} && $params_to_app->{ref_type} ne "none") {
        open(my $afa, '<', "$work_dir/$prefix.afa") or die "Cannot open $work_dir/$prefix.afa: $!";
        my $ref_seq = "";
        my $ref_count = 0;
        while ( my $line = <$afa> ) {
                chomp($line); # remove newlines
                $line =~ s/#.*//; # remove comments
                $line =~ s/;.*//; # remove comments
                $line =~ s/^\s+//;  # remove leading whitespace
                $line =~ s/\s+$//; # remove trailing whitespace
                next if(length($line) <= 0);
                if(substr($line, 0, 1) eq ">" && $ref_count > 0) {
                    last;
                } elsif (substr($line, 0, 1) eq ">") {
                    $ref_count += 1;
                    next;
                }
                $ref_seq = $ref_seq . $line;
        }
        close($afa);
        open(my $table, '<', "$work_dir/foma.table") or die "Cannot open $work_dir/foma.table: $!";
        open(my $new_table, '>', "$work_dir/foma2.table") or die "Cannot open $work_dir/foma2.table: $!";
        my $first_line = 1;
        my $pos_count = 0;
        my $spaces_count = 0;
        while ( my $line = <$table> ) {
            if ($first_line == 1) {
                $first_line = 0;
                print $new_table snp_line_mod($line, "Position in Reference", "Reference Character");
            } else {
                my $char = substr($ref_seq, $pos_count, 1);
                $pos_count += 1;
                if ($char ne "-") {
                    my $curr_pos = $pos_count - $spaces_count;
                    print $new_table snp_line_mod($line, $curr_pos, $char);
                } else {
                    print $new_table snp_line_mod($line, "", $char);
                    $spaces_count += 1;
                }
            }
        }
        close($table);
        close($new_table);
        rename "$work_dir/foma2.table", "$work_dir/$prefix.snp.tsv";
        unlink("$work_dir/foma.table") or warn "Unable to unlink: $work_dir/foma.table: $!";
    } else {
        rename "$work_dir/foma.table", "$work_dir/$prefix.snp.tsv";
    }
    #
    # Make a gene tree.
    #
    my $out_type = "aligned_protein_fasta";
    my $tree_alphabet = "protein";
    if ($dna) {
        $out_type = "aligned_dna_fasta";
        $tree_alphabet = "DNA";
    }
    @cmd = ("p3x-build-gene-tree", "--program", "fasttree", "--alphabet", $tree_alphabet, "--output_dir", "$work_dir", "$work_dir/$prefix.afa");
    run_cmd(\@cmd);
    print STDOUT "Completed gene tree creation.\n";

    # 
    # Reroot the fasttree nwk file using midpoint rooting, then delete the fasttree nwk
    #
    my $fasttree_nwk = "$work_dir/${prefix}_fasttree.nwk";
    # make sure initial nwk file exists
    if (-e $fasttree_nwk) {
        eval {
            my $midpoint_tmp = "$work_dir/${prefix}.nwk";
            my $midpoint_nwk = "$work_dir/${prefix}_midpoint.nwk";
            my @midpoint_cmd = ("p3x-reformat-tree","--midpoint","-f","newick","-o",$midpoint_tmp,"-i",$fasttree_nwk);
            print "midpoint root command: @midpoint_cmd\n";
            run_cmd(\@midpoint_cmd);
            unlink $fasttree_nwk;
            # rename midpoint file to same name as fasttree file
            # ensures MSA tree viewer does not break
            # terminate if rename doesn't work?:or die "Unable to rename: $!"
            rename($midpoint_nwk,$fasttree_nwk);
        };
        if ($@) {
            print "error while midpoint rooting nwk file: \n$@\n";
        }
    }

    #
    # Copy output to the workspace.
    #
    my @output_suffixes = (
        [qr/\.afa$/, $out_type],
        [qr/\.nexus$/, "txt"],
        [qr/\.phy$/, "txt"],
        [qr/\.pir$/, "txt"],
        [qr/\.xmfa$/, "txt"],
        [qr/\.mauve\.log$/, "txt"],
        [qr/_log\.txt$/, "txt"],
        [qr/\.nwk$/, "nwk"],
        [qr/\.job.log$/, "txt"],
        [qr/\.aln$/, "txt"],
        [qr/\.consensus\.fasta$/, "txt"],
        [qr/\.tsv$/, "tsv"],
        # [qr/\.table$/, "tsv"],
        [qr/\.png$/, "png"],
        [qr/\.svg$/, "svg"],
        );
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    my @files = sort { $a cmp $b } grep { -f "$work_dir/$_" } readdir(D);
    my $output = 1;
    for my $file (@files)
    {
	for my $suf (@output_suffixes)
	{
	    if ($file =~ $suf->[0])
	    {
 	    	$output=0;
		my $path = "$output_folder/$file";
		my $type = $suf->[1];
		$app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_folder/$file", $type, 1,
					       (-s "$work_dir/$file" > 10_000 ? 1 : 0), # use shock for larger files
					       $token);
	    }
	}
    }
    #
    # Clean up staged input files.
    #
    while (my($orig, $staged_file) = each %$staged)
    {
	    unlink($staged_file) or warn "Unable to unlink $staged_file: $!";
    }
    unlink($text_input_file) or warn "Unable to unlink $text_input_file: $!";
    unlink($ofile) or warn "Unable to unlink $ofile: $!";
    return $output;
}

sub run_cmd() {
    my $cmd = $_[0];
    my $ok = run(@$cmd);
    if (!$ok)
    {
        die "Command failed: @$cmd\n";
    }
}

sub is_aa {
    my ($file_str, $is_file) = @_;
    $is_file //= 1;
    if ($is_file) {
        open(FH, "<", $file_str) or die "Could not open file $file_str to check if it has amino acids in it. $!";
    } else {
        open(FH, "<", \$file_str) or die "Could not open string $file_str to check if it has amino acids in it. $!";
    }
    my $dna_count = 0;
    my $str_len = 0;
    while (my $line = <FH>) {
        if ((substr($line, 0, 1) ne ">")) { # and not($line =~ /^[ACTGNactgn-]+$/)
            $line =~ tr/-//;
            my $str_len += length($line);
            my $dna_count += $line =~ tr/ACTGNactgn//;
        } else {
            if ($str_len > 0 && ($dna_count / $str_len < $THRESHOLD)) {
                close FH or die $!;
                return 1;
            }
            $dna_count = 0;
            $str_len = 0;
        }
    }
    if ($str_len > 0 && ($dna_count / $str_len < $THRESHOLD)) {
        close FH or die $!;
        return 1;
    }
    close FH or die $!;
    return 0;
}

sub convert_aa_line {
    my ($line) = @_;
    chomp($line);
    $line = uc $line;
    my @codons = unpack '(A3)*', $line;
    my @aminoAcids = map { exists $aacode{$_} ? $aacode{$_} : "X" } @codons;
    return join('', @aminoAcids);
}

sub convert_aa_file {
    my($in_file, $out_file, $is_file) = @_;
    print STDERR "Converting a file to amino acids.\n";
    if ($is_file) {
        open(INF, "<", $in_file) or die "Couldn't open file $in_file. $!";
    } else {
        open(INF, "<", \$in_file) or die "Couldn't open string. $!";
    }
    open(OUTF, ">", $out_file) or die "Couldn't open file $out_file. $!";
    while (my $line = <INF>) {
        if (substr($line, 0, 1) eq ">") {
            print OUTF "$line";
        } else {
            print OUTF  convert_aa_line($line) . "\n";
        }
    }
    close INF or die $!;
    close OUTF or die $!;
}

sub get_genome_group_file {
    my($data_api_module, $genome_group, $target_dir, $exclude_id) = @_;
    my $ids = $data_api_module->retrieve_patric_ids_from_genome_group($genome_group);
    my $filename = basename($genome_group);
    my $work_fasta = "$target_dir/$filename.fasta";
    open(my $in, '>', $work_fasta) or die "Cannot open $work_fasta: $!";
    get_genome_seqs($data_api_module, $ids, $in, $target_dir, $exclude_id);
    close($in);
    return $work_fasta;
}

sub get_genome_seqs {
    my($data_api_module, $ids, $in, $target_dir, $exclude_id) = @_;
    $data_api_module->retrieve_contigs_in_genomes($ids, $target_dir, "%s");
    for my $gid (@$ids) {
        next if $gid eq $exclude_id;
        my $loc = "$target_dir/$gid";
        open my $fh, '<', $loc or die "Cannot open $loc: $!";
        my $seq_line = "";
        my $count_contigs = 0;
        while ( my $line = <$fh> ) {
            # chomp($line); # remove newlines
            # $line =~ s/#.*//; # remove comments
            # $line =~ s/;.*//; # remove comments
            # $line =~ s/^\s+//;  # remove leading whitespace
            # $line =~ s/\s+$//; # remove trailing whitespace
            next if(length($line) <= 0);
            if (substr($line, 0, 1) ne ">") {
                $seq_line = $seq_line . $line;
            } elsif ($count_contigs > 0) {
                last;
            } else {
                $count_contigs = $count_contigs + 1;
                print $in ">" . $gid . " " . substr($line, 1);
            }
        }
        print $in $seq_line;
        close($loc);
        unlink($loc) or warn "Unable to unlink $loc: $!";
    }
}

sub snp_line_mod {
    my($line, $first, $second) = @_;
    chomp($line);
    my @spl = split("\t", $line);
    my $line_str = "";
    my $arr_cnt = 0;
    for my $i (@spl) {
        $line_str = $line_str . "\t" . $i;
        $arr_cnt += 1;
        if ($arr_cnt == 1) {
            $line_str = $line_str . "\t" . $first;
        } elsif ($arr_cnt == 2) {
            $line_str = $line_str . "\t" . $second;
        }
    }
    $line_str =~ s/^\s+//;
    return $line_str . "\n";
}
