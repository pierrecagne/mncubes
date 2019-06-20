#!/usr/bin/perl -W
# Compute ./mncubes m n mu with fixed m n for certain values of mu
# and return the total number of m,n-cubes.

use POSIX;

foreach $opt (@ARGV) {
    if($opt eq "-h") {
	print
"$0 : help
This program counts the number of (m,n)-cubes given m and n.
It requires the program ./mncubes to work.

Usage -- $0 [-h] m n

Options --
  -h : will display this help and quit (other options or 
       arguments are ignored)
";
	exit 0;
    }
}

($m, $n) = @ARGV;
die "Wrong usage. Use $0 -h for more informations.\n" unless defined($m) && defined($n);

$den_max = ($m-1)*($n-1);
($gen_val) = `./mncubes $m $n 1/$den_max` =~ /M_.*=\s*(\d+)/;

# $href_farey = &table_farey($den_max-1);
# $sum = $href_farey->{"0/1"};

# for (keys %$href_farey) {
#     $sum += $gen_val - $href_farey->{$_};
# }

# print "M_{$m,$n} with Farey = $sum \n";

print "mu = 1/$den_max : $gen_val \n";

$sum = 0;
for (my $i = 2; $i < $den_max; ++$i) {
    $sum += &totient($i);
}
$sum +=1 ; $sum *= $gen_val;
for (my $i = 2; $i < $den_max; ++$i) {
    my ($tmp) = `./mncubes $m $n 1/$i` =~ /M_.*=\s*(\d+)/;
    print "mu = 1/$i : $tmp \n";
    $sum -= &totient($i)*$tmp;
}
print "M_{$m,$n} = $sum \n";

sub totient {
    my ($k) = @_;
    my $res = 0;
    for (my $i = 1; $i < $k; ++$i) {
	$res++ if &gcd($k,$i) == 1;
    }
    return $res;
}

# sub table_farey {
#     my ($max) = @_;
#     my ($a,$b,$c,$d) = (0,1,1,$max);
#     my $href_table = {};
    
#     while ($c < $max) {
# 	if($a <= 1) {
# 	    ($href_table->{"$a/$b"}) = `./mncubes $m $n $a/$b` =~ /M_.*=\s*(\d+)/;
# 	} else {
# 	    $href_table->{"$a/$b"} = $href_table->{"1/$b"};
# 	}
# 	print "mu = $a/$b : ".$href_table->{"$a/$b"}." \n";
# 	# my ($num,$den) = ($a*$d+$b*$c,2*$b*$d);
# 	# $gcd = &gcd ($num,$den);
# 	# $num /= $gcd; $den /= $gcd;
# 	# my ($tmp) = `./mncubes $m $n $num/$den` =~ /M_.*=\s*(\d+)/;
# 	# print "mu = $num/$den : $tmp \n"; 
# 	$k = floor(($max+$b)/$d);
# 	($a,$b,$c,$d) = ($c, $d, $k*$c-$a, $k*$d-$b);
#     }

#     return $href_table;
# }

sub gcd {
    my ($a,$b) = @_;
    &gcd ($b, $a) unless $a >= $b;
    if($b == 0) {
	return $a;
    } else {
	&gcd ($b, $a % $b);
    }
}
