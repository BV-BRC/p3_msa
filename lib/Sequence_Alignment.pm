package Sequence_Alignment;

use strict;
use warnings;

use List::Util ('any');
use IPC::Run 'run';

our $debug = 0;

sub new {
    my ($class, $input) = @_;
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {};
    $self->{_ids} = [];
    $self->{_loci} = [];
    $self->{_default_locus} = '';
    $self->{_length} = {};
    $self->{_metadata} = {};

    if ($input) {
        $self->read_file($input);
    }
    return $self;
}

sub set_debug {
    my $onoff = shift;
    $debug = $onoff ? 1 : 0;
}

sub get_ntaxa {
    my $self = shift;
    return scalar(@{$self->{_ids}});
}

sub get_length {
    my ($self, $locus) = @_;
    if ($locus) {
        return $self->{_length}{$locus};
    }

    my $total_length = 0;
    for my $loc (@{$self->{_loci}}) {
        $total_length += $self->{_length}{$loc};
    }
    return $total_length;
}

sub get_ids {
    my ($self, $locus) = @_;
    if ($locus) {
        return keys(%{$self->{_seqs}{$locus}});
    }
    return @{$self->{_ids}};
}

sub get_locus_ids {
    my $self = shift;
    return @{$self->{_loci}};
}

sub get_sequence_lengths {
    my $self = shift;
    my %lengths;
    for my $locus (@{$self->{_loci}}) {
        for my $id (keys %{$self->{_seqs}{$locus}}) {
            $lengths{$id}{$locus} = length($self->{_seqs}{$locus}{$id});
        }
    }
    return \%lengths;
}

sub add_seq {
    my ($self, $id, $seq, $locus) = @_;
    if ($locus) {
        if (!$self->{_default_locus}) {
            $self->{_default_locus} = $locus;
        }
    }
    else {
        if (!$self->{_default_locus}) {
            $self->{_default_locus} = 'default';
        }
        $locus = $self->{_default_locus};
    }

    if (!exists $self->{_seqs}{$locus}) {
        $self->{_seqs}{$locus} = {};
        $self->{_length}{$locus} = 0;
        push @{$self->{_loci}}, $locus;
        my @sorted_loci = sort(@{$self->{_loci}});
        $self->{_loci} = \@sorted_loci;
    }

    if (!$self->{_in_alignment} && exists $self->{_seqs}{$locus}{$id}) {
        warn("duplicate occurrence of id $id at locus $locus");
        my $temp = $id;
        my $suffix = 1;
        while (exists $self->{_seqs}{$locus}{$temp}) {
            $suffix++;
            $temp = "${id}_$suffix";
        }
        $id = $temp;
    }

    $seq =~ tr/ //d;
    $self->{_seqs}{$locus}{$id} = $seq;
    if (length($seq) > $self->{_length}{$locus}) {
        $self->{_length}{$locus} = length($seq);
    }
    unless (grep { $_ eq $id } @{$self->{_ids}}) {
        push @{$self->{_ids}}, $id;
    }
    return $id;
}

sub set_metadata {
    my ($self, $id, $category, $value) = @_;
    if (!exists $self->{_metadata}{$category}) {
        $self->{_metadata}{$category} = {};
    }
    $self->{_metadata}{$category}{$id} = $value;
}

sub get_metadata {
    my ($self, $id, $category) = @_;
    if (exists $self->{_metadata}{$category}) {
        return $self->{_metadata}{$category}{$id};
    }
    return undef;
}

sub detect_file_format {
    my ($self, $fh) = @_;
    $_ = readline $fh;
    seek $fh, 0, 0;

    my $format = 'unknown';
    $format = 'clustal' if (/^CLUSTAL/ || /^MUSCLE/);
    $format = 'fasta' if (/^>/);
    $format = 'phylip' if (/^(\d+)\s+(\d+)\s*$/);
    $format = 'nexus' if (/\#NEXUS/);
    return $format;
}

sub read_file {
    my ($self, $fh, $format, $locus) = @_;
    my @retval;
    my $filename = $fh;
    my $opened_here = 0;

    if (!ref($fh)) {
        open(my $opened_fh, '<', $filename) or die "Cannot open $filename: $!";
        $fh = $opened_fh;
        $opened_here = 1;
        if (!$locus) {
            $locus = $filename;
            $locus =~ s/.*\///;
            $locus =~ s/\..*//;
        }
    }
    if (!$locus) {
        $locus = $self->{_default_locus};
    }
    if (!$format) {
        $format = $self->detect_file_format($fh);
    }
    if ($format eq 'unknown') {
        close($fh) if $opened_here;
        return undef;
    }

    if ($format eq 'clustal') {
        my $found = 0;
        while (<$fh>) {
            if (/^CLUSTAL/ || /^MUSCLE/) {
                $found = 1;
                last;
            }
        }
        die "Format seems to be wrong, not Clustal.\n" if !$found;
        while (<$fh>) {
            if (/^(\S+)\s+(\S+)/) {
                my ($id, $seq) = ($1, $2);
                $id = $self->add_seq($id, $seq, $locus);
                push @retval, $id;
            }
        }
    }
    elsif ($format eq 'phylip') {
        $_ = <$fh>;
        die "Format does not seem to be phylip\n" unless /^\s*(\d+)\s+(\d+)\s*$/;
        my $ntaxa = $1;
        my %seqhash;
        my @ids;
        for my $i (1 .. $ntaxa) {
            $_ = <$fh>;
            /(\S+)\s+(\S.*\S)/ or die $_;
            my ($id, $seq) = ($1, $2);
            push @ids, $id;
            $seqhash{$id} = $seq;
        }
        my $index = 0;
        while (<$fh>) {
            my $seq = $_;
            $seq =~ s/\s//g;
            if ($seq) {
                my $id = $ids[$index % $ntaxa];
                $seqhash{$id} .= $seq;
                $index += 1;
            }
        }
        for my $id (@ids) {
            $seqhash{$id} =~ s/\./\-/g;
            $id = $self->add_seq($id, $seqhash{$id}, $locus);
            push @retval, $id;
        }
    }
    elsif ($format eq 'fasta') {
        my ($id, $seq) = ('', '');
        while (<$fh>) {
            if (/^>(\S+)/) {
                if ($seq) {
                    $id = $self->add_seq($id, $seq, $locus);
                    push @retval, $id;
                    $seq = '';
                }
                $id = $1;
            }
            elsif (/(.*\S)/) {
                $seq .= $1;
            }
        }
        if ($seq) {
            $id = $self->add_seq($id, $seq, $locus);
            push @retval, $id;
        }
    }
    elsif ($format eq 'nexus') {
        $_ = <$fh>;
        die "Format does not seem to be NEXUS" unless /\#NEXUS/;
        my $matrix;
        while (<$fh>) {
            chomp;
            s/\[[^\]]*\]//g;
            if (/^matrix/i) {
                $matrix = 1;
                next;
            }
            if ($matrix) {
                if (/^(\S+)\s+(\S+)/) {
                    my ($id, $seq) = ($1, $2);
                    $id = $self->add_seq($id, $seq, $locus);
                    push @retval, $id;
                }
                last if /;/;
            }
        }
    }
    close($fh) if $opened_here;
    return @retval;
}

sub write_to_file {
    my ($self, $fh, $format, $locus) = @_;
    my $out;
    if (ref($fh)) {
        $out = $fh;
    }
    else {
        open($out, '>', $fh) or die "Cannot open $fh: $!";
    }
    $format ||= 'fasta';
    my @ids = $self->get_ids();
    my %seqs;
    my @loci;

    if (defined($locus) && any { $_ eq $locus } @{$self->{_loci}}) {
        @loci = ($locus);
    }
    else {
        @loci = $self->get_locus_ids();
    }

    for my $loc (@loci) {
        for my $id (@ids) {
            if (exists($self->{_seqs}{$loc}{$id}) && $self->{_seqs}{$loc}{$id}) {
                $seqs{$id} .= $self->{_seqs}{$loc}{$id};
            }
            else {
                $seqs{$id} .= '-' x $self->get_length($loc);
            }
        }
    }

    my $max_id_length = 0;
    if ($format eq 'phylip') {
        print $out $self->get_ntaxa(), "  ", $self->get_length(), "\n";
        for my $id (@ids) {
            if (length($id) > $max_id_length) {
                $max_id_length = length($id);
            }
        }
    }

    for my $id (@ids) {
        if (($format eq 'fasta') or ($format eq 'raxml')) {
            if ($self->{_in_alignment}) {
                $seqs{$id} =~ tr/-//d;
            }
            my $output_id = $id;
            if ($format eq 'raxml') {
                my $changed = $output_id =~ tr/:()[]/_____/;
                if ($changed) {
                    if (!exists($self->{_raxml_to_original_id})) {
                        $self->{_raxml_to_original_id} = {};
                    }
                    $self->{_raxml_to_original_id}{$output_id} = $id;
                }
            }
            print $out ">$output_id\n$seqs{$id}\n";
        }
        elsif ($format eq 'phylip') {
            printf($out "%-${max_id_length}s  %s\n", $id, $seqs{$id});
        }
    }

    if (!ref($fh)) {
        close $out;
    }
}

sub write_fasta {
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'fasta', $locus);
}

sub write_phylip {
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'phylip', $locus);
}

sub write_fasta_for_raxml {
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'raxml', $locus);
}

sub restore_original_ids_in_raxml_tree {
    my ($self, $newick) = @_;
    return $newick unless exists $self->{_raxml_to_original_id};
    for my $raxml_id (keys %{$self->{_raxml_to_original_id}}) {
        $newick =~ s/$raxml_id/$self->{_raxml_to_original_id}{$raxml_id}/;
    }
    return $newick;
}

sub align {
    my ($self, $locus, $options) = @_;
    my @command_list;
    my $details = "";
    $options ||= {};

    if ($locus) {
        my $aligner = lc($options->{aligner} || 'mafft');
        my $strategy = lc($options->{strategy} || 'auto');
        my $work_dir = $options->{work_dir} || '.';
        my $unaligned = "$work_dir/${locus}_unaligned.fa";
        my $aligned_fasta = "$work_dir/${locus}_aligned.fa";
        my $stderr = "$work_dir/${locus}_alignment_stderr.txt";
        unlink($stderr);

        $self->write_to_file($unaligned, 'fasta', $locus);

        my $cmd;
        if ($aligner eq 'muscle') {
            $cmd = ['muscle', '-quiet', '-in', $unaligned, '-fastaout', $aligned_fasta];
        }
        elsif ($aligner eq 'mafft') {
            my $parallel = $ENV{P3_ALLOCATED_CPU};
            if ($strategy eq 'fftns1') {
                $cmd = ['mafft', '--retree', 1, '--maxiterate', 0];
            }
            elsif ($strategy eq 'fftns2') {
                $cmd = ['fftns'];
            }
            elsif ($strategy eq 'fftnsi') {
                $cmd = ['fftnsi'];
            }
            elsif ($strategy eq 'einsi') {
                $cmd = ['einsi'];
            }
            elsif ($strategy eq 'linsi') {
                $cmd = ['linsi'];
            }
            elsif ($strategy eq 'ginsi') {
                $cmd = ['ginsi'];
            }
            else {
                $cmd = ['mafft', '--auto', '--preservecase'];
            }
            if ($parallel && ($cmd->[0] eq 'mafft')) {
                push @$cmd, '--thread', $parallel;
            }
            push @$cmd, $unaligned;
        }
        else {
            die "Unsupported aligner '$aligner' for locus-aware sequence alignment.\n";
        }

        my $retval;
        if ($aligner eq 'muscle') {
            $retval = run($cmd, "2>", $stderr);
        }
        else {
            $retval = run($cmd, ">", $aligned_fasta, "2>", $stderr);
        }
        die "Alignment command failed for locus $locus: @$cmd\n" unless $retval;

        $self->{_in_alignment} = 1;
        $self->read_file($aligned_fasta, 'fasta', $locus);
        delete($self->{_in_alignment});

        if (open(my $err, '<', $stderr)) {
            local $/;
            $details .= join(" ", @$cmd) . "\n" . <$err>;
            close $err;
        }
        else {
            $details .= join(" ", @$cmd) . "\n";
        }
        push @command_list, join(" ", @$cmd) . "\n";
    }
    else {
        for my $loc (@{$self->{_loci}}) {
            my ($cmd_aref, $detail) = $self->align($loc, $options);
            push @command_list, @$cmd_aref;
            $details .= $detail;
        }
    }

    return \@command_list, $details;
}

sub calc_column_gap_count {
    my ($self, $locus) = @_;
    my @gap_count;
    $#gap_count = $self->{_length}{$locus} - 1;
    for my $id (keys %{$self->{_seqs}{$locus}}) {
        my @str_as_array = split('', $self->{_seqs}{$locus}{$id});
        for my $gap_pos (0 .. $#str_as_array) {
            $gap_count[$gap_pos] += $str_as_array[$gap_pos] eq '-';
        }
    }
    return \@gap_count;
}

sub end_trim {
    my ($self, $threshold, $locus) = @_;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my @steps;
    my $details = "";

    if ($locus) {
        my $gap_count = $self->calc_column_gap_count($locus);
        my $max_gaps = (1.0 - $threshold) * scalar(keys %{$self->{_seqs}{$locus}});
        my $left_trim = 0;
        my $right_trim = 0;
        $left_trim++ while (($gap_count->[$left_trim] > $max_gaps) and $left_trim < $self->{_length}{$locus} - 1);
        $right_trim++ while (($gap_count->[$self->{_length}{$locus} - 1 - $right_trim] > $max_gaps) and $right_trim < $self->{_length}{$locus} - 1);
        $right_trim++ while (($right_trim > $left_trim) and ($gap_count->[$right_trim - 1] > $max_gaps));
        if ($left_trim or $right_trim) {
            my $len = $self->{_length}{$locus} - $right_trim - $left_trim;
            push @steps, "Trim $left_trim on left, trim $right_trim on right, new length = $len.\n";
            for my $id (keys %{$self->{_seqs}{$locus}}) {
                $self->{_seqs}{$locus}{$id} = substr($self->{_seqs}{$locus}{$id}, $left_trim, $len);
            }
            $self->{_length}{$locus} = $len;
        }
        else {
            push @steps, "No trimming needed for $locus\n";
        }
        $details = join("", @steps);
    }
    else {
        for my $loc (@{$self->{_loci}}) {
            my ($locus_steps, $locus_detail) = $self->end_trim($threshold, $loc);
            push @steps, @{$locus_steps};
            $details .= $locus_detail;
        }
    }
    return \@steps, $details;
}

sub calc_row_gap_count {
    my ($self, $locus) = @_;
    my %gap_count;
    for my $id (keys %{$self->{_seqs}{$locus}}) {
        $gap_count{$id} = $self->{_seqs}{$locus}{$id} =~ tr/-/-/;
    }
    return \%gap_count;
}

sub get_alignment_occupancy {
    my $self = shift;
    my %occupancy;
    for my $locus (@{$self->{_loci}}) {
        my $length = $self->{_length}{$locus};
        my $gaps_per_seq = $self->calc_row_gap_count($locus);
        my $num_seqs = scalar keys %{$gaps_per_seq};
        my $total_gaps = 0;
        $total_gaps += $_ for values %{$gaps_per_seq};
        my $total_chars = $length * $num_seqs;
        my $percent_non_gap = $total_chars ? 100 * (1 - $total_gaps / $total_chars) : 0;
        $occupancy{$locus} = { length => $length, percent_non_gap => $percent_non_gap };
    }
    return \%occupancy;
}

sub delete_gappy_seqs {
    my ($self, $threshold, $locus) = @_;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my @retval;

    if ($locus) {
        my $max_gaps = (1.0 - $threshold) * $self->get_length($locus);
        for my $id (keys %{$self->{_seqs}{$locus}}) {
            my $num_gaps = $self->{_seqs}{$locus}{$id} =~ tr/-/-/;
            if ($num_gaps > $max_gaps) {
                delete($self->{_seqs}{$locus}{$id});
                my $command = "sequence $id";
                $command .= " at segment $locus" unless $locus eq 'default';
                $command .= sprintf(" deleted due to $num_gaps gaps (%.2f%%), exceeding threshold of %.2f%%\n", $num_gaps / $self->get_length($locus), $threshold);
                push @retval, $command;
            }
        }
        if (scalar @retval) {
            $self->{_length}{$locus} = 0;
            for my $id (keys %{$self->{_seqs}{$locus}}) {
                if (length($self->{_seqs}{$locus}{$id}) > $self->{_length}{$locus}) {
                    $self->{_length}{$locus} = length($self->{_seqs}{$locus}{$id});
                }
            }
        }
    }
    else {
        for my $loc (@{$self->{_loci}}) {
            my $deleted_seqs = $self->delete_gappy_seqs($threshold, $loc);
            if ($deleted_seqs && scalar(@{$deleted_seqs})) {
                push @retval, @{$deleted_seqs};
            }
        }
        unless (scalar @retval) {
            @retval = ("No sequences gappy enough to remove");
        }
    }
    return \@retval;
}

sub write_stats {
    my $self = shift;
    my %prop_gaps;
    for my $locus (@{$self->{_loci}}) {
        my $gaps_per_seq = $self->calc_row_gap_count($locus);
        my $total_gaps = 0;
        for my $id (keys %{$gaps_per_seq}) {
            $total_gaps += $gaps_per_seq->{$id};
        }
        $prop_gaps{$locus} = $total_gaps / ((scalar keys %{$gaps_per_seq}) * $self->{_length}{$locus});
    }

    my $retval = "";
    $retval .= "Alignment Statistics   " . join("\t", @{$self->{_loci}}) . "\n";
    $retval .= "Number of sequences: ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . scalar(keys %{$self->{_seqs}{$locus}});
    }
    $retval .= "\n";

    $retval .= "Alignment length:    ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . $self->get_length($locus);
    }
    $retval .= "\n";

    $retval .= "Proportion gaps:     ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . sprintf("%.4f", $prop_gaps{$locus});
    }
    $retval .= "\n";

    return $retval;
}

1;
