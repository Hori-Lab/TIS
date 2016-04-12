#!/usr/bin/perl -w

$define = "" ;
$ext    = "F90" ;

# make directory "cpp/"
$dir = "./cpp/" ;
unless(-d $dir) {
   mkdir($dir,0755) ;
}

# C Pre Processor
@file = glob "*.$ext" ;
foreach (@file)
{
   system("cpp $define $_ > $dir$_") ;
}

exit ;
