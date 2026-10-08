#The MSA application with variance analysis.

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use feature "switch";
use FindBin;
use lib "$FindBin::Bin/../lib";
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
use SeedUtils;
use gjoseqlib qw();
use Path::Tiny qw(path);
use Sequence_Alignment;

use experimental qw(switch);


my $script = Bio::KBase::AppService::AppScript->new(\&process_fasta, \&preflight);
my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;

# Threshold for determining if input data is AA or DNA
my $THRESHOLD = 0.75;

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight 
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $token = $app->token();
    my $ws = $app->workspace();

    my $api = P3DataAPI->new();

    # get number of genomes from genome groups, if any
    my $numGenomes = 0;
    my $genomegroups = defined($params->{select_genomegroup}) ? $params->{select_genomegroup} : undef;
    if (defined $genomegroups) 
    {
        for my $gg (@$genomegroups)
        {
            print "$gg\n";
            my $genomes = retrieve_genome_group_ids($api, $ws, $gg);
            my $n = @$genomes;
            $numGenomes = $numGenomes + $n;
        }
    }
    if (exists($params->{genome_list})) {
        $numGenomes = $numGenomes + scalar(@{$params->{genome_list}})
    }

    # get number of features from feature groups
    my $numFeatures = 0;
    my $featuregroups = defined($params->{feature_groups}) ? $params->{feature_groups} : undef;
    if (defined $featuregroups) 
    {
        for my $fg (@$featuregroups)
        {
            print "$fg\n";
            my $features;
            if ($params->{alphabet} eq "protein")
            {
                $features = $api->retrieve_protein_sequences_from_feature_group($fg);
            } 
            else 
            {
                $features = $api->retrieve_nucleotide_sequences_from_feature_group($fg);
            }
            if (defined $features)
            {
                my $nf = @$features;
                $numFeatures = $numFeatures + $nf;
            }
        }
    }
    if (exists($params->{feature_list})) {
        $numFeatures = $numFeatures + scalar(@{$params->{feature_list}})
    }

    # get file sizes
    my $files = defined($params->{fasta_files}) ? $params->{fasta_files} : undef;
    my $totFileSize = 0;
    if (defined $files) 
    {
        for my $ff (@$files)
        {
            my $path = $ff->{file};
            if (!defined($path))
            {
                die "Preflight: a fasta_files entry is missing its 'file' path: " . Dumper($ff);
            }
            print "$path\n";
            my $file_data = $ws->stat($path);
            if (!defined($file_data))
            {
                die "Preflight: cannot stat input fasta file '$path' " .
                    "(it does not exist, or the submitting user lacks permission to read it)\n";
            }
            $totFileSize = $totFileSize + $file_data->size;
        }
    } 

    my $input_type = defined($params->{input_type}) ? $params->{input_type} : undef;

    my $runtime = 0;
    my $mem = '';
    my $mem_threshold = 50000000000; #50GB
    if (defined $input_type and $input_type eq "input_sequence")
    {
        $runtime = "1800";
        $mem = "8GB";
    }
    else 
    {
        my $numGroups = $numFeatures + $numGenomes;
        if ($numGroups == 0) 
        {
            $runtime = 3 * 3600;    
            $mem = '32GB';
        } elsif ($numGroups < 10 and $totFileSize < ($mem_threshold/10)) {
            $runtime = 1800; 
            $mem = '8GB';
        } elsif ($numGroups < 100 and $totFileSize < ($mem_threshold/5))  {
            $runtime = 3 * 3600;
            $mem = '16GB';
        } elsif ($numGroups < 500 and $totFileSize < ($mem_threshold/2)) {
            $runtime = 6 * 3600;
            $mem = '32GB';
        } elsif ($numGroups < 1000 and $totFileSize < ($mem_threshold/2)) {
            $runtime = 43200; # 12 hours
            $mem = '32GB';
        } elsif ($numGroups < 3000 and $totFileSize >= $mem_threshold) {
            $runtime = 43200 * 2;
            $mem = '64GB';
        } else { # <= 5000 genomes
            $runtime = 43200 * 4;
            $mem = '64GB';
        } 
    }
    # zero genome ids: default
    # have no reference for this so just guessing
    my $pf = {
        cpu => 1,
        memory => $mem,
        runtime => $runtime,
        storage => 0,
        is_control_task => 0
    };
    return $pf;
}

sub process_fasta
{
    my($app, $app_def, $raw_params, $params) = @_;
    print "Proc MSA Var ", Dumper($app_def, $raw_params, $params);
    my $token_object = $app->token();
    my $token = ref($token_object) && $token_object->can('token')
        ? $token_object->token()
        : $token_object;
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
    if (exists($params_to_app->{feature_list})) {
        $file_count = $file_count + 1;
    }
    if (exists($params_to_app->{genome_list})) {
        $file_count = $file_count + 1;
    }
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
    my $genome_selection = $params_to_app->{genome_selection} && ref($params_to_app->{genome_selection}) eq 'HASH'
        ? $params_to_app->{genome_selection}
        : {};
    my $selected_segments = $genome_selection->{selected_segments} && ref($genome_selection->{selected_segments}) eq 'ARRAY'
        ? $genome_selection->{selected_segments}
        : [];
    my $concat_segments = $genome_selection->{concat_segments} ? 1 : 0;
    my $segment_aware_genome_groups = (scalar(@$selected_segments) > 0)
        && (exists($params_to_app->{select_genomegroup}) || exists($params_to_app->{genome_list}));
    my $has_input_fasta_files = scalar(@{$params_to_app->{fasta_files} || []}) > 0;
    my $has_feature_groups = scalar(@{$params_to_app->{feature_groups} || []}) > 0;
    my $has_feature_list = scalar(@{$params_to_app->{feature_list} || []}) > 0;
    my $has_keyboard_input = defined($params_to_app->{fasta_keyboard_input})
        && length($params_to_app->{fasta_keyboard_input}) >= 1;
    if (exists($params_to_app->{select_genomegroup})) {
        $dna = 1;
        $in_type = "feature_dna_fasta";
        $file_count = $file_count + scalar(@{$params_to_app->{select_genomegroup}});
        unless ($segment_aware_genome_groups) {
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
    }
    say STDERR "Number of input sources: $file_count.";
    if (exists($params_to_app->{genome_list})) {
        $dna = 1;
        $in_type = "feature_dna_fasta";
        unless ($segment_aware_genome_groups) {
            my $genome_list = $params_to_app->{genome_list};
            my $genome_id_exclude = "";
            if ($params_to_app->{ref_type} eq "genome_id") {
                $genome_id_exclude = $params_to_app->{ref_string};
                $genome_id_exclude =~ s/^\s+|\s+$//g;
            }
            push @genome_groups, get_genome_list_file($data_api_module, $genome_list, $stage_dir, $genome_id_exclude);
        }
    }
    #
    # Write files to the staging directory.
    #
    my $aligned_exists = 0;
    my $mixed = 0;
    my $staged = {};
    my $ofile = "$stage_dir/feature_groups.fasta";
    my $text_input_file = "$stage_dir/fasta_keyboard_input.fasta";
    unless ($segment_aware_genome_groups) {
        my @to_stage;
        for my $read_tuple (@{$params_to_app->{fasta_files} || []}) {
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
        open(F, ">$ofile") or die "Could not open $ofile";
        my $feature_id_exclude = "";
        if ($params_to_app->{ref_type} eq "feature_id") {
            $feature_id_exclude = $params_to_app->{ref_string};
            $feature_id_exclude =~ s/^\s+|\s+$//g;
        }
        for my $feature_name (@{$params_to_app->{feature_groups} || []}) {
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
        if ($has_feature_list) {
            my $feature_list = $params_to_app->{feature_list};
            my @ids_new = ();
            for my $id (@$feature_list){
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
        if ($has_feature_groups || $has_feature_list) {
	    push @{ $params_to_app->{fasta_files} }, {"file" => $ofile, "type" => $in_type};
        }
        close(F);
        #
        # Put keyboard input into a file.
        #
        if ($has_keyboard_input) {
            if ((not $dna) && (not is_aa($params_to_app->{fasta_keyboard_input}, 0))) {
                convert_aa_file($params_to_app->{fasta_keyboard_input}, $text_input_file, 0);
            } else {
                open(FH, '>', $text_input_file) or die "Cannot open $text_input_file: $!";
                print FH $params_to_app->{fasta_keyboard_input};
                close(FH);
            }
            push @{ $params_to_app->{fasta_files} }, {"file" => $text_input_file, "type" => $in_type};
        }
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
            for my $read_tuple (@{$params_to_app->{fasta_files} || []}) {
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
    }
    my $recipe = lc($params_to_app->{aligner});
    my @rand_set = ('0' ..'9', 'A' .. 'F');
    my %id_map = ();
    my @analysis_sets;
    my $primary_alignment_written = 0;
    my $segment_genome_count = 0;
    if ($segment_aware_genome_groups) {
        if ($has_input_fasta_files || $has_feature_groups || $has_feature_list || $has_keyboard_input) {
            die "Genome segment selection for MSA currently supports genome group or genome list inputs only.\n";
        }
        if ($recipe eq 'progressivemauve' or $recipe eq 'mauve') {
            die "Segment-aware genome group alignment currently supports Muscle and Mafft only.\n";
        }
        my $seq_al = build_segmented_genome_group_alignment($data_api_module, $app->workspace(), $params_to_app, $selected_segments, $stage_dir);
        my @loci = $seq_al->get_locus_ids();
        my $genome_count = $seq_al->get_ntaxa();
        $segment_genome_count = $genome_count;
        die "No genome sequences matched the selected segments: " . join(", ", @$selected_segments) . ".\n"
            unless @loci && $genome_count;
        die "At least two genomes are required for segment-aware MSA; found $genome_count.\n"
            if $genome_count < 2;
        say STDERR "Segment-aware input: $genome_count genomes across " . scalar(@loci) . " selected segment loci.";
        my %locus_sequence_counts;
        for my $locus (@loci) {
            my @ids = $seq_al->get_ids($locus);
            $locus_sequence_counts{$locus} = scalar(@ids);
        }
        $seq_al->align(undef, {
            aligner => $recipe,
            strategy => lc($params_to_app->{strategy} // 'auto'),
            work_dir => $work_dir
        });
        if ($params_to_app->{trim_threshold}) {
            $seq_al->end_trim($params_to_app->{trim_threshold});
        }
        if ($params_to_app->{gap_threshold}) {
            $seq_al->delete_gappy_seqs($params_to_app->{gap_threshold});
        }
        for my $locus (@loci) {
            my $segment_folder = $locus;
            $segment_folder =~ s/[^A-Za-z0-9_.-]+/_/g;
            my $segment_dir = "$work_dir/$segment_folder";
            -d $segment_dir or mkdir $segment_dir or die "Cannot mkdir $segment_dir: $!";

            my $segment_base = "${prefix}_${segment_folder}";
            $seq_al->write_fasta("$segment_dir/${segment_base}.afa", $locus);
            push @analysis_sets, {
                base => $segment_base,
                local_dir => $segment_dir,
                result_folder => $segment_folder,
                sequence_count => $locus_sequence_counts{$locus},
                expected_sequence_count => $genome_count,
                alignment_length => $seq_al->get_length($locus)
            };
        }
        if ($concat_segments) {
            $seq_al->write_fasta("$work_dir/output.afa");
            $primary_alignment_written = 1;
            unshift @analysis_sets, {
                base => $prefix,
                local_dir => $work_dir,
                result_folder => "",
                sequence_count => $genome_count,
                expected_sequence_count => $genome_count,
                alignment_length => $seq_al->get_length()
            };
        }
    }
    else {
        #
        # Combine all files into one input.fasta file.
        # Put the reference sequence first if present.
        #
        # TODO: download feature id sequences and put them into input_fasta
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
            print STDERR "CHECK $read_tuple->{type} $filename convert=$convert\n";
            open my $fh, '<', $filename or die "Cannot open $filename: $!";
            my $seq_line = "";


            print STDERR "Translating $filename to $work_fasta\n";
            while (my($id, $def, $seq) = gjoseqlib::read_next_fasta_seq($fh))
            {
                my $print_me = 1;

                my $unique_id = join '' => map $rand_set[rand @rand_set], 1 .. 8;
                $id_map{$id} = $unique_id;
                print STDERR "$id is replaced by $unique_id\n";

                # Remove indels from alignments if other files are present.
                if ($aligned_exists && $file_count > 1)
                {
                    $seq =~ tr/-_.~*//d;
                }
                if ($convert)
                {
                    my $aa = translate(\$seq);
                    gjoseqlib::print_alignment_as_fasta(\*IN, [$unique_id, $def, $aa]);
                }
                else
                {
                    gjoseqlib::print_alignment_as_fasta(\*IN, [$unique_id, $def, $seq]);
                }
            }
            close($fh);
        }
        close(IN);
        #
        # Run the multiple sequence aligner.
        #
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
        my $strategy = lc($params_to_app->{strategy});
        print STDOUT "Selected strategy: $strategy";
        my @mafft_cmd;
        given ($strategy) {
            when ($_ eq "fftns1") {
                @mafft_cmd = ("mafft", "--retree", 1, "--maxiterate", 0, "$work_dir/input.fasta");
            }
            when ($_ eq "fftns2") {
                @mafft_cmd = ("fftns", "$work_dir/input.fasta");
            }
            when ($_ eq "fftnsi") {
                @mafft_cmd = ("fftnsi", "$work_dir/input.fasta");
            }
            when ($_ eq "einsi") {
                @mafft_cmd = ("einsi", "$work_dir/input.fasta");
            }
            when ($_ eq "linsi") {
                @mafft_cmd = ("linsi", "$work_dir/input.fasta");
            }
            when ($_ eq "ginsi") {
                @mafft_cmd = ("ginsi", "$work_dir/input.fasta");
            }
            default {
                @mafft_cmd = ("mafft", "--auto", "--preservecase", "$work_dir/input.fasta");
            }	
        } 
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
        $primary_alignment_written = 1;
        push @analysis_sets, {
            base => $prefix,
            local_dir => $work_dir,
            result_folder => ""
        };
    }
    if ($primary_alignment_written) {
        rename "$work_dir/output.afa", "$work_dir/$prefix.afa"
            or die "Unable to rename output alignment to $prefix.afa: $!";
    }
    my $out_type = "aligned_protein_fasta";
    my $tree_alphabet = "protein";
    if ($dna) {
        $out_type = "aligned_dna_fasta";
        $tree_alphabet = "DNA";
    }
    say STDERR "Alignment result sets: " . join(", ", map { $_->{base} } @analysis_sets) . ".";
    for my $analysis_set (@analysis_sets) {
        process_alignment_outputs($analysis_set->{local_dir}, $analysis_set->{base}, $dna, $tree_alphabet, $params_to_app);
    }
    write_msa_report(
        "$work_dir/index.html",
        $output_folder,
        $prefix,
        \@analysis_sets,
        $params,
        $selected_segments,
        $concat_segments,
        $segment_aware_genome_groups,
        $segment_genome_count,
        $dna
    );

    #
    # Copy output to the workspace.
    #
    my @output_suffixes = (
        [qr/\.afa$/, $out_type, 1],
        [qr/\.nexus$/, "txt", 1],
        [qr/\.phy$/, "txt", 1],
        [qr/\.pir$/, "txt", 1],
        [qr/\.xmfa$/, "txt", 0],
        [qr/\.mauve\.log$/, "txt", 0],
        [qr/_log\.txt$/, "txt", 0],
        [qr/\.nwk$/, "nwk", 1],
        [qr/\.job.log$/, "txt", 0],
        [qr/\.aln$/, "txt", 1],
        [qr/\.consensus\.fasta$/, "txt", 0],
        [qr/\.tsv$/, "tsv", 0],
        # [qr/\.table$/, "tsv", 0],
        [qr/\.png$/, "png", 0],
        [qr/\.svg$/, "svg", 0],
        [qr/\.html$/, "html", 0],
        );
    my @upload_locations = ({
        local_dir => $work_dir,
        workspace_dir => $output_folder,
        result_folder => ""
    });
    my %seen_result_folders;
    for my $analysis_set (@analysis_sets) {
        my $result_folder = $analysis_set->{result_folder};
        next unless $result_folder && !$seen_result_folders{$result_folder}++;

        my $workspace_dir = "$output_folder/$result_folder";
        my $created = eval {
            $app->workspace->create({
                overwrite => 1,
                objects => [[$workspace_dir, "folder"]]
            });
            1;
        };
        die "Unable to create workspace result folder $workspace_dir: $@" unless $created;
        push @upload_locations, {
            local_dir => $analysis_set->{local_dir},
            workspace_dir => $workspace_dir,
            result_folder => $result_folder
        };
    }

    my $output = 1;
    for my $upload_location (@upload_locations) {
        my $local_dir = $upload_location->{local_dir};
        my $workspace_dir = $upload_location->{workspace_dir};
        opendir(my $output_dir, $local_dir) or die "Cannot opendir $local_dir: $!";
        my @files = sort { $a cmp $b } grep { -f "$local_dir/$_" } readdir($output_dir);
        closedir($output_dir);

        for my $file (@files) {
	    for my $suf (@output_suffixes) {
	        if ($file =~ $suf->[0]) {
		    $output=0;
		    my $type = $suf->[1];
                my $is_replace_id = $suf->[2];
                my $local_file = "$local_dir/$file";
                if ($is_replace_id) {
                    restore_id($local_file, %id_map);
                }
                my $saved = eval {
                    $app->workspace->save_file_to_file($local_file, {}, "$workspace_dir/$file", $type, 1,
                                                       (-s $local_file > 10_000 ? 1 : 0), # use shock for larger files
                                                       $token);
                    1;
                };
                if ($saved) {
                    say STDERR "Saved $file to $workspace_dir as $type.";
                }
                else {
                    warn "Workspace upload failed for $file; retrying with p3-cp: $@";
                    my ($extension) = $file =~ /\.([^.]+)$/;
                    $extension ||= "dat";
                    my @copy_cmd = ("p3-cp", "-f", "-m", "$extension=$type", $local_file, "ws:$workspace_dir");
                    run_cmd(\@copy_cmd);
                    say STDERR "Saved $file to $workspace_dir as $type using p3-cp.";
                }
	        }
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
    if ($has_keyboard_input && -e $text_input_file) {
        unlink($text_input_file) or warn "Unable to unlink $text_input_file: $!";
    }
    if (($has_feature_groups || $has_feature_list) && -e $ofile) {
        unlink($ofile) or warn "Unable to unlink $ofile: $!";
    }
    return $output;
}

sub write_msa_report {
    my ($report_file, $workspace_dir, $prefix, $analysis_sets, $params, $selected_segments,
        $concat_segments, $segment_aware, $genome_count, $dna) = @_;

    open(my $report, '>', $report_file) or die "Cannot write MSA report $report_file: $!";

    my @input_sources;
    push @input_sources, @{ $params->{select_genomegroup} || [] };
    push @input_sources, @{ $params->{feature_groups} || [] };
    push @input_sources, map { $_->{file} || () } @{ $params->{fasta_files} || [] };
    push @input_sources, scalar(@{ $params->{genome_list} || [] }) . ' genome IDs'
        if $params->{genome_list} && @{ $params->{genome_list} };
    push @input_sources, scalar(@{ $params->{feature_list} || [] }) . ' feature IDs'
        if $params->{feature_list} && @{ $params->{feature_list} };
    push @input_sources, 'Pasted sequence data'
        if defined($params->{fasta_keyboard_input}) && length($params->{fasta_keyboard_input});

    my $inputs = @input_sources
        ? join('<br>', map { html_escape($_) } @input_sources)
        : 'Not available';
    my $segments = @$selected_segments
        ? join(', ', map { html_escape($_) } @$selected_segments)
        : 'Not applicable';
    my $concat = $segment_aware ? ($concat_segments ? 'Yes' : 'No') : 'Not applicable';
    my $aligner = html_escape($params->{aligner} || 'Not specified');
    my $strategy = html_escape($params->{strategy} || 'auto');
    my $trim_threshold = defined($params->{trim_threshold}) ? html_escape($params->{trim_threshold}) : 'Not specified';
    my $gap_threshold = defined($params->{gap_threshold}) ? html_escape($params->{gap_threshold}) : 'Not specified';
    my $alphabet = $dna ? 'DNA' : 'Protein';

    print $report <<"HTML";
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Multiple Sequence Alignment Report</title>
  <style>
    body { color: #333; font: 14px Arial, sans-serif; line-height: 1.45; margin: 24px; }
    h1, h2 { color: #24557a; }
    h1 { margin-bottom: 6px; }
    .subtitle { color: #666; margin: 0 0 24px; }
    table { border-collapse: collapse; margin-bottom: 24px; width: 100%; }
    th, td { border: 1px solid #c9c9c9; padding: 8px 10px; text-align: left; vertical-align: top; }
    th { background: #edf2f5; }
    .summary th { width: 190px; }
    .warning { background: #fff6df; border-left: 4px solid #c58a00; margin-bottom: 20px; padding: 10px 12px; }
    .muted { color: #777; }
    a { color: #075985; }
  </style>
</head>
<body>
  <h1>Multiple Sequence Alignment Report</h1>
  <p class="subtitle">Output: @{[html_escape($prefix)]}</p>
  <h2>Run Summary</h2>
  <table class="summary">
    <tr><th>Input</th><td>$inputs</td></tr>
    <tr><th>Alphabet</th><td>$alphabet</td></tr>
    <tr><th>Alignment program</th><td>$aligner</td></tr>
    <tr><th>Strategy</th><td>$strategy</td></tr>
    <tr><th>Selected segments</th><td>$segments</td></tr>
    <tr><th>Concatenate segments</th><td>$concat</td></tr>
    <tr><th>Trim threshold</th><td>$trim_threshold</td></tr>
    <tr><th>Gap threshold</th><td>$gap_threshold</td></tr>
  </table>
HTML

    if ($segment_aware && !$concat_segments) {
        print $report '<div class="warning">No concatenated alignment was requested. Results are organized by segment.</div>';
    }

    print $report <<"HTML";
  <h2>Results</h2>
  <table>
    <thead>
      <tr>
        <th>Analysis</th>
        <th>Sequences present</th>
        <th>Missing genomes</th>
        <th>Alignment length</th>
        <th>Files</th>
      </tr>
    </thead>
    <tbody>
HTML

    for my $analysis (@$analysis_sets) {
        my $folder = $analysis->{result_folder} || '';
        my $base = $analysis->{base};
        my $label = $folder ? $folder : ($segment_aware ? 'Combined' : 'Alignment');
        $label =~ s/^segment_/Segment /;
        my $sequence_count = defined($analysis->{sequence_count})
            ? $analysis->{sequence_count}
            : fasta_sequence_count("$analysis->{local_dir}/$base.afa");
        my $expected = $analysis->{expected_sequence_count} || $genome_count || $sequence_count;
        my $missing = $expected > $sequence_count ? $expected - $sequence_count : 0;
        my $alignment_length = defined($analysis->{alignment_length})
            ? $analysis->{alignment_length}
            : fasta_alignment_length("$analysis->{local_dir}/$base.afa");
        my $relative_base = $folder ? "$folder/$base" : $base;
        my @links;
        for my $file_spec (
            ['Alignment', '.afa'],
            ['Tree', '_fasttree.nwk'],
            ['SNP table', '.snp.tsv'],
            ['Consensus', '.consensus.fasta'],
            ['Entropy plot', '.entropy.svg']
        ) {
            my ($file_label, $suffix) = @$file_spec;
            my $local_file = "$analysis->{local_dir}/$base$suffix";
            next unless -f $local_file;
            push @links, workspace_file_link($workspace_dir, "$relative_base$suffix", $file_label);
        }

        print $report '<tr><td>' . html_escape($label) . '</td>';
        print $report '<td>' . html_escape($sequence_count) . '</td>';
        print $report '<td>' . html_escape($missing) . '</td>';
        print $report '<td>' . html_escape($alignment_length) . '</td>';
        print $report '<td>' . (@links ? join(' | ', @links) : '<span class="muted">No files found</span>') . '</td></tr>';
    }

    print $report <<"HTML";
    </tbody>
  </table>
  <p class="muted">All generated files are also available in the workspace result folder.</p>
</body>
</html>
HTML
    close($report) or die "Cannot close MSA report $report_file: $!";
}

sub workspace_file_link {
    my ($workspace_dir, $relative_path, $label) = @_;
    my $path = "$workspace_dir/$relative_path";
    $path =~ s{//+}{/}g;
    my $href = '/workspace' . join('/', map { uri_escape($_) } split('/', $path, -1));
    return '<a href="' . html_escape($href) . '" target="_blank">' . html_escape($label) . '</a>';
}

sub html_escape {
    my ($value) = @_;
    $value = '' unless defined($value);
    $value =~ s/&/&amp;/g;
    $value =~ s/</&lt;/g;
    $value =~ s/>/&gt;/g;
    $value =~ s/"/&quot;/g;
    $value =~ s/'/&#39;/g;
    return $value;
}

sub fasta_sequence_count {
    my ($file) = @_;
    return 0 unless -f $file;
    open(my $fh, '<', $file) or return 0;
    my $count = 0;
    while (my $line = <$fh>) {
        $count++ if $line =~ /^>/;
    }
    close($fh);
    return $count;
}

sub fasta_alignment_length {
    my ($file) = @_;
    return 0 unless -f $file;
    open(my $fh, '<', $file) or return 0;
    my $length = 0;
    my $current = 0;
    while (my $line = <$fh>) {
        if ($line =~ /^>/) {
            $length = $current if $current > $length;
            $current = 0;
        }
        else {
            $line =~ s/\s+//g;
            $current += length($line);
        }
    }
    $length = $current if $current > $length;
    close($fh);
    return $length;
}

sub process_alignment_outputs {
    my ($work_dir, $analysis_base, $dna, $tree_alphabet, $params_to_app) = @_;
    my $alignment_file = "$work_dir/${analysis_base}.afa";
    say STDERR "Processing result set for $analysis_base.";

    open(my $alignment_fh, '<', $alignment_file) or die "Cannot open $alignment_file: $!";
    open(my $output_fh, '>', "$work_dir/output.afa") or die "Cannot open $work_dir/output.afa: $!";
    while (my $line = <$alignment_fh>) {
        next unless length($line);
        print $output_fh substr($line, 0, 1) eq '>' ? $line : uc($line);
    }
    close($alignment_fh);
    close($output_fh);

    my @cmd = ("snp_analysis", "-r", $work_dir, "-x");
    push @cmd, "-n" if $dna;
    run_cmd(\@cmd);
    unlink("$work_dir/output.afa") or warn "Unable to unlink $work_dir/output.afa: $!";
    print STDOUT "Completed SNP analysis for $analysis_base.\n";
    rename "$work_dir/cons.fasta", "$work_dir/${analysis_base}.consensus.fasta"
        or warn "Unable to save consensus for $analysis_base: $!";

    @cmd = ("convert_seq_files", $alignment_file, "$work_dir/$analysis_base", $dna ? "DNA" : "protein");
    run_cmd(\@cmd);

    @cmd = ("snp_analysis_figure", "$work_dir/foma.table", "$work_dir/$analysis_base");
    run_cmd(\@cmd);
    print STDOUT "Completed figure creation for $analysis_base.\n";

    if ($params_to_app->{ref_type} && $params_to_app->{ref_type} ne "none") {
        open(my $afa, '<', $alignment_file) or die "Cannot open $alignment_file: $!";
        my $ref_seq = "";
        my $ref_count = 0;
        while (my $line = <$afa>) {
            chomp($line);
            $line =~ s/#.*//;
            $line =~ s/;.*//;
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            next unless length($line);
            if (substr($line, 0, 1) eq ">" && $ref_count > 0) {
                last;
            }
            elsif (substr($line, 0, 1) eq ">") {
                $ref_count++;
                next;
            }
            $ref_seq .= $line;
        }
        close($afa);

        open(my $table, '<', "$work_dir/foma.table") or die "Cannot open $work_dir/foma.table: $!";
        open(my $new_table, '>', "$work_dir/foma2.table") or die "Cannot open $work_dir/foma2.table: $!";
        my $first_line = 1;
        my $pos_count = 0;
        my $spaces_count = 0;
        while (my $line = <$table>) {
            if ($first_line) {
                $first_line = 0;
                print $new_table snp_line_mod($line, "Position in Reference", "Reference Character");
                next;
            }
            my $char = substr($ref_seq, $pos_count, 1);
            $pos_count++;
            if ($char ne "-") {
                my $curr_pos = $pos_count - $spaces_count;
                print $new_table snp_line_mod($line, $curr_pos, $char);
            }
            else {
                print $new_table snp_line_mod($line, "", $char);
                $spaces_count++;
            }
        }
        close($table);
        close($new_table);
        rename "$work_dir/foma2.table", "$work_dir/${analysis_base}.snp.tsv"
            or die "Unable to save SNP table for $analysis_base: $!";
        unlink("$work_dir/foma.table") or warn "Unable to unlink $work_dir/foma.table: $!";
    }
    else {
        rename "$work_dir/foma.table", "$work_dir/${analysis_base}.snp.tsv"
            or die "Unable to save SNP table for $analysis_base: $!";
    }

    @cmd = ("p3x-build-gene-tree", "--program", "fasttree", "--alphabet", $tree_alphabet,
            "--output_dir", ".", basename($alignment_file));
    run_cmd_in_dir(\@cmd, $work_dir);
    print STDOUT "Completed gene tree creation for $analysis_base.\n";

    my $fasttree_nwk = "$work_dir/${analysis_base}_fasttree.nwk";
    if (-e $fasttree_nwk) {
        eval {
            my $midpoint_nwk = "$work_dir/${analysis_base}_midpoint.nwk";
            my @midpoint_cmd = ("p3x-reformat-tree", "--midpoint", "-f", "newick",
                                "-o", $midpoint_nwk, "-i", $fasttree_nwk);
            print "midpoint root command: @midpoint_cmd\n";
            run_cmd(\@midpoint_cmd);
            unlink($fasttree_nwk) or die "Unable to remove unrooted tree $fasttree_nwk: $!";
            rename($midpoint_nwk, $fasttree_nwk)
                or die "Unable to save midpoint-rooted tree $fasttree_nwk: $!";
        };
        warn "Error while midpoint rooting $analysis_base: $@" if $@;
    }
}

sub run_cmd {
    my $cmd = $_[0];
    my $ok = run(@$cmd);
    if (!$ok)
    {
        die "Command failed: @$cmd\n";
    }
}

sub run_cmd_in_dir {
    my ($cmd, $dir) = @_;
    my $original_dir = Cwd::getcwd();
    chdir($dir) or die "Cannot change directory to $dir: $!";

    my $ok = eval {
        run_cmd($cmd);
        1;
    };
    my $error = $@;

    chdir($original_dir) or die "Cannot restore working directory to $original_dir: $!";
    die $error unless $ok;
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
    while (my($id, $def, $seq) = gjoseqlib::read_next_fasta_seq(\*FH))
    {
	$seq =~ tr/-//;
	$str_len += length($seq);
	$dna_count += $seq =~ tr/ACTGNactgn//;
    }
    close(FH);
    return $str_len > 0 && ($dna_count / $str_len < $THRESHOLD);
}

sub restore_id {
  my($filepath, %id_map) = @_;

  my $file = path($filepath);
  my $data = $file->slurp_utf8;

  while( my($id, $unique_id) = each %id_map)
  {
    $data =~ s/$unique_id/$id/g;
  }
  $file->spew_utf8($data);
}

sub convert_aa_file {
    my($in_file, $out_file, $is_file) = @_;
    if ($is_file) {
	print STDERR "Convert $in_file to AA sequences in $out_file\n";
        open(INF, "<", $in_file) or die "Couldn't open file $in_file. $!";
    } else {
	my $frag = substr($in_file, 0, 10);
	$frag =~ s/\r//g;
	print STDERR "Convert DNA starting with $frag to AA sequences in $out_file\n";
        open(INF, "<", \$in_file) or die "Couldn't open string. $!";
    }
    open(OUTF, ">", $out_file) or die "Couldn't open file $out_file. $!";

    while (my($id, $def, $seq) = gjoseqlib::read_next_fasta_seq(\*INF))
    {
	my $aa = translate(\$seq);
	gjoseqlib::print_alignment_as_fasta(\*OUTF, [$id, $def, $aa]);
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

sub get_genome_list_file {
    my($data_api_module, $genome_list, $target_dir, $exclude_id) = @_;
    my $work_fasta = "$target_dir/genome_list.fasta";
    open(my $in, '>', $work_fasta) or die "Cannot open $work_fasta: $!";
    get_genome_seqs($data_api_module, $genome_list, $in, $target_dir, $exclude_id);
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

sub build_segmented_genome_group_alignment {
    my ($data_api_module, $workspace, $params_to_app, $selected_segments, $target_dir) = @_;
    my $seq_al = Sequence_Alignment->new();
    my %selected = map { $_ => 1 } @{$selected_segments || []};

    my $genome_id_exclude = "";
    if ($params_to_app->{ref_type} && $params_to_app->{ref_type} eq "genome_id") {
        $genome_id_exclude = $params_to_app->{ref_string};
        $genome_id_exclude =~ s/^\s+|\s+$//g;
    }

    if ($genome_id_exclude) {
        add_segmented_genome_sequences($seq_al, $data_api_module, [$genome_id_exclude], \%selected, "reference", "");
    }
    if (exists($params_to_app->{select_genomegroup})) {
        for my $genome_group (@{$params_to_app->{select_genomegroup}}) {
            my $genome_ids = retrieve_genome_group_ids($data_api_module, $workspace, $genome_group);
            die "Genome group $genome_group did not contain any genome IDs.\n"
                unless $genome_ids && @$genome_ids;
            say STDERR "Segment-aware genome group $genome_group contains " . scalar(@$genome_ids) . " genomes.";
            add_segmented_genome_sequences($seq_al, $data_api_module, $genome_ids, \%selected, basename($genome_group), $genome_id_exclude);
        }
    }
    if (exists($params_to_app->{genome_list})) {
        add_segmented_genome_sequences($seq_al, $data_api_module, $params_to_app->{genome_list}, \%selected, "genome_list", $genome_id_exclude);
    }

    return $seq_al;
}

sub retrieve_genome_group_ids {
    my ($data_api_module, $workspace, $genome_group) = @_;
    my $genome_ids = $data_api_module->retrieve_patric_ids_from_genome_group($genome_group);
    return $genome_ids if $genome_ids && @$genome_ids;

    my $objects = eval {
        $workspace->get({
            objects => [$genome_group],
            metadata_only => 0
        });
    };
    if ($@) {
        warn "Unable to read genome group $genome_group from the workspace: $@";
        return [];
    }
    return [] unless $objects && ref($objects) eq 'ARRAY' && @$objects;

    my $data = $objects->[0]->[1];
    if (defined($data) && !ref($data)) {
        $data = eval { decode_json($data) };
        if ($@) {
            warn "Unable to parse genome group $genome_group data: $@";
            return [];
        }
    }

    return [] unless $data
        && ref($data) eq 'HASH'
        && $data->{id_list}
        && ref($data->{id_list}) eq 'HASH'
        && $data->{id_list}{genome_id}
        && ref($data->{id_list}{genome_id}) eq 'ARRAY';

    my @ids = grep { defined($_) && $_ ne '' } @{$data->{id_list}{genome_id}};
    say STDERR "Read " . scalar(@ids) . " genome IDs directly from workspace group $genome_group.";
    return \@ids;
}

sub add_segmented_genome_sequences {
    my ($seq_al, $data_api_module, $genome_ids, $selected, $source_label, $exclude_id) = @_;

    for my $genome_id (@{$genome_ids || []}) {
        next if $exclude_id && $genome_id eq $exclude_id;

        my ($resp, $data) = $data_api_module->submit_query('genome_sequence', "eq(genome_id,$genome_id)", "sequence,contig,segment");
        die "Genome sequence query failed for genome_id $genome_id.\n"
            unless $data && ref($data) eq 'ARRAY';
        say STDERR "Genome $genome_id returned " . scalar(@$data) . " sequence records.";
        my %segments_seen;
        my $matched_records = 0;
        for my $record (@$data) {
            next unless $record->{sequence};

            my $segment = defined($record->{segment}) ? $record->{segment} : '';
            $segments_seen{$segment ne '' ? $segment : '(none)'} = 1;
            my $use = !%{$selected || {}};
            if (%{$selected || {}} && $segment ne '') {
                $use = exists $selected->{$segment};
            }
            next unless $use;

            my $locus = $segment ne '' ? "segment_$segment" : undef;
            my $id = $seq_al->add_seq($genome_id, $record->{sequence}, $locus);
            $seq_al->set_metadata($id, 'data_source', $source_label);
            $seq_al->set_metadata($id, 'database_link', 'genome_id');
            $matched_records++;
        }
        say STDERR "Genome $genome_id matched $matched_records records; segments seen: "
            . join(", ", sort keys %segments_seen) . ".";
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
