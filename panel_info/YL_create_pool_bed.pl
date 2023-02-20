use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long 'GetOptions';
use List::Util qw[min max sum];


sub grepPoolBed {
	my $bedFile = shift;
	my $poolId  = shift;
	my $outputBedGrep = shift;
	my $outputBedFull = shift;

	my @amplicon;
	open (BEDFILE, "<$bedFile") or fatality("Can't open BED file");
	open (WRITEBED, ">$outputBedFull");
	
	## Цикл ниже открывает входной bed файл, заданных опцией '-t' и ищет все ампликоны и добавляет их в переменную @amplicon
	## Структура элемента массива @amplcion: [хромосома, позиция начала, позиция конца, номер пула]
	while (<BEDFILE>) {
		chomp;
		my $line = $_;
		my $str = $line;
		my @mas = split/\t/;
		$str = uc($str);
		if ($str =~ /POOL=(\d+)/) {
			my $pool = $1;
			if ((defined($amplicon[0]))
				and($amplicon[scalar @amplicon - 1]->[0] eq $mas[0])
				and($amplicon[scalar @amplicon - 1]->[1] > $mas[1])) {
				die "Input panel targets Bed file is not Sorted; use sort -k1,1 -k2,2n\n";
				}
			push(@amplicon, [$mas[0], $mas[1], $mas[2], $pool]);
			if ($pool eq $poolId) {
				print WRITEBED "$line\n";
				}
			} else {
			die "Wrong Input Bed File Format. Please include Pool=XX into each string to defined amplicon pool\n";
			}
		}
	close WRITEBED;
	close BEDFILE;
	open (BEDWRITE, ">$outputBedGrep");
	
	## Цикл ниже бежит по всем элементам из массива @amplicon
	## Для каждого ампликона из пула PoolId (входная опция сабрутины) определяется точка наиболее отсоящая от всех позиций начала и конца пула
	for (my $link = 0; $link < scalar @amplicon; $link++) {
		next unless $amplicon[$link]->[3] eq $poolId;
		my $distance = 0;
		my $position;
		for (my $l = $amplicon[$link]->[1]; $l <= $amplicon[$link]->[2]; $l++) {
			my $min = 9999999;
			for (my $link_i = $link - 10; $link_i < $link + 10; $link_i++){
				next unless defined $amplicon[$link_i];
				next if $amplicon[$link_i]->[0] ne $amplicon[$link]->[0];
				my $d1 = abs($amplicon[$link_i]->[1] - $l);
				my $d2 = abs($amplicon[$link_i]->[2] - $l);
				my $d = min($d1, $d2);
				if ($d < $min) {$min = $d}
				}
			if ($min > $distance) {
				$distance = $min;
				$position = $l;
				}
			}
		print BEDWRITE "$amplicon[$link]->[0]\t",$position-1,"\t",$position,"\n";
		}
	close BEDWRITE;
	}

sub base_name {
	my $name = shift;
	if ($name =~ /(\S+).bed$/) {
		return $1;
		}
	return undef;
	}

sub run {
	my ($options) = @_;
	open (my $bed, "<", $options->{bed});
	
	print STDERR base_name($options->{bed}),"\n";
	my %pool_DIC;
	while (<$bed>) {
		chomp;
		my $str = $_;
		$str = uc($str);
		if ($str =~ /POOL=(\d+)/) {
			$pool_DIC{$1} = 1;
			} else {
			die "Wrong Input Bed File Format. Please include Pool=XX into each string to defined amplicon pool\n";
			}
		}
	close $bed;
	
	foreach my $pool (keys %pool_DIC) {
		grepPoolBed($options->{bed}, $pool, 
			base_name($options->{bed}) . ".pool_$pool.grep.bed",
			base_name($options->{bed}) . ".pool_$pool.full.bed");
		}
	}

sub option_builder {
        my ($factory) = @_;
        my %opts;

        &GetOptions (
                'h|help'        => \$opts{'h'},
                'b|bed=s'       => \$opts{'bed'},
        );
        pod2usage(0) if($opts{'h'});
        pod2usage(1) if(!$opts{'bed'});

        return \%opts;
}


{
  my $options = option_builder();
  run($options);
}


__END__

=head1 NAME

Grep pool bed

=head1 SYNOPSIS

This program creates new bed file based on the given with pool identifying positions.

    -b   input bed file

=cut

