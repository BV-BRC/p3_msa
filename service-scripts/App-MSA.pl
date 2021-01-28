#The MSA application with variance analysis.

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use lib "/homes/jsporter/p3_msa/p3_core/lib";
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
    if (exists $params_to_app{fasta_files}) {
        $file_count = $file_count + length($params_to_app{fasta_files});
    }
    if (exists $params_to_app{feature_groups}) {
        $file_count = $file_count + length($params_to_app{feature_groups});
    }
    if (length($params_to_app->{fasta_keyboard_input}) >= 1) {
        $file_count = $file_count + 1;
    }
    say STDOUT "Number of files: $file_count.";
    my $prefix = $params_to_app->{output_file};
    #
    # Determine if the data is represented as DNA or protein.
    #
    my $dna = 1;
    my $in_type = "feature_dna_fasta";
    if (substr($params_to_app->{alphabet}, 0, 1) eq "p") {
    	$dna = 0;
	$in_type = "feature_protein_fasta"
    }
    #
    # Write files to the staging directory.
    #
    my @to_stage;
    my $aligned_exists = 0;
    for my $read_tuple (@{$params_to_app->{fasta_files}})
    {
	for my $read_name (keys %{$read_tuple})
	{
	   if($read_name eq "file")
           {
	       my $nameref = \$read_tuple->{$read_name};
	       $in_files{$$nameref} = $nameref;
	       push(@to_stage, $$nameref);
           }
        if(rindex $read_name, "aligned", 0) {
            $aligned_exists = 1;
        }
        }
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
    for my $feature_name (@{$params_to_app->{feature_groups}}) {
	    my $ids = $data_api_module->retrieve_patricids_from_feature_group($feature_name);
	    my $seq = "";
	    if ($dna) {
		$seq = $data_api_module->retrieve_nucleotide_feature_sequence($ids);
	    } else {
		$seq = $data_api_module->retrieve_protein_feature_sequence($ids);
	    }
	    for my $id (@$ids) {
		    my $out = ">$id\n" . $seq->{$id} . "\n";
    		    print F $out;
	    }
    }
    if (exists($params_to_app->{feature_groups})) {
    	push @{ $params_to_app->{fasta_files} }, {"file" => $ofile, "type" => $in_type};
	close(F);
	# delete $params_to_app->{feature_groups};
    }
    #
    # Put keyboard input into a file.
    #
    my $text_input_file = "$stage_dir/fasta_keyboard_input.fasta";
    open(FH, '>', $text_input_file) or die "Cannot open $text_input_file: $!";
    print FH $params_to_app->{fasta_keyboard_input};
    push @{ $params_to_app->{fasta_files} }, {"file" => $text_input_file, "type" => $in_type};
    close(FH);
    # Combine all files into one input.fasta file.
    my $work_fasta = "$work_dir/input.fasta";
    open(IN, '>', $work_fasta) or die "Cannot open $work_fasta: $!";
    for my $read_tuple (@{$params_to_app->{fasta_files}}) {
    	my $filename = $read_tuple->{file};
	open my $fh, '<', $filename or die "Cannot open $filename: $!";
	while ( my $line = <$fh> ) {
	        chomp; # remove newlines
            if ($aligned_exists && $file_count > 1 && not substr($line, 0, 1) ne ">") {
                s/[-._*]+//g;
            }
                s/#.*//; # remove comments
                s/;.*//; # remove comments
                s/^\s+//;  # remove leading whitespace
                s/\s+$//; # remove trailing whitespace
                next if(length($line) <= 0);
		print IN $line;
	}
	close($fh);
    }
    close(IN);
    #
    # Run the multiple sequence aligner.
    #
    if ($aligned_exists && $file_count == 1) {
        rename "$work_dir/input.fasta", "$work_dir/$prefix.afa";
    }
    elsif (lc($params_to_app{aligner}) eq "muscle") {
    	my $muscle_cmd =  "muscle -in $work_dir/input.fasta -fastaout $work_dir/$prefix.afa -clwout $work_dir/$prefix.aln";
    } else {
        die "Recipe not found: $recipe\n";
    }
    # Run the SNP analysis.
    my @cmd = ("/homes/jsporter/p3_msa/p3_msa/lib/snp_analysis.pl", "-r", "$work_dir");
    if ($dna) {
    	push @cmd, "-n";
    }
    run_cmd(\@cmd);
    rename "$work_dir/cons.fasta", "$work_dir/$prefix.cons.fasta";
    #
    # Create figures.
    #
    @cmd = ("python3", "/homes/jsporter/p3_msa/p3_msa/lib/snp_analysis_figure.py", "$work_dir/foma.table", "$work_dir/$prefix");
    run_cmd(\@cmd);
    #
    # Copy output to the workspace.
    #
    rename "$work_dir/foma.table", "$work_dir/$prefix.tsv";
    my $out_type = "aligned_protein_fasta";
    if ($dna) {
        $out_type = "aligned_dna_fasta";
    }
    my @output_suffixes = (
        [qr/\.afa$/, $out_type],
        [qr/\.aln$/, "txt"],
        [qr/\cons.fasta$/, "txt"],
        [qr/\.tsv$/, "tsv"],
        [qr/\.table$/, "tsv"],
        [qr/\.png$/, "png"],
        [qr/\.svg$/, "svg"]);
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    my @files = sort { $a cmp $b } grep { -f "$work_dir/$_" } readdir(D);
    my $output=1;
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
