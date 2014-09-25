#!/usr/bin/perl
$filename = "phi-base-1_vs36.txt";
open (MYFILE, $filename) || die "Sorry, could open $filename\n";
while (<MYFILE>) {
   ($a,$b,$c,$d,$e,$f,$g,$h,$i,$j,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z,$aa,$ab,$ac,$ad,$ae,$af,$ag,$ah,$ai,$aj,$ak,$al,$am,$an,$ao,$ap,$aq,$ar,$as,$at,$au,$av,$aw,$ax) = split(/\t/,$_);
   write; # write to stdout
}
close (MYFILE);

format STDOUT =
@<<<<<<<<< @<<<<<<<<<<<<<<< @<<<<<<<<<<<< @<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<<<<<<<<<<<
$a,        $b,              $c,           $d,            $e,                 $f,                   $g
.
