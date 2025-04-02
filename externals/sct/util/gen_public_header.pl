#!/usr/bin/perl

use strict;
use warnings;

my @h=();

sub ph {
  my @x = @_;
  for my $t (@x) {
    chomp $t;
    push @h, "$t\n";
  }
}

sub proc_file {
  my $f = shift;
  open(F,"<$f") || die "cannot open $f";
  my $public = 0;
  my $x;
  ph("// from $f:");
  while($x=<F>) {
    if ( $x =~ /^\s*\/\/\s*public([+-])\s*$/ ) {
      my $s = $1;
      if ($s eq "+") {
        die "+/- mismatch" if $public;
        $public=1;
        next;
      } elsif ($s eq "-") {
        die "+/- mismatch" unless $public;
        $public=0;
        next;
      }
      print "x=[$x] => public=$public\n";
    }
    ph($x) if $public;
  }
  close(F);
}


if (!@ARGV) {
  print STDERR "usage: $0 <header files>\n";
  exit(1);
}

ph("/* sct.h: do not edit - edit $0 instead */");
ph("#ifndef SCT_HEADER_FILE_USED");
ph("#define SCT_HEADER_FILE_USED");
my $mpidef = shift @ARGV;
ph("#ifndef $mpidef");
ph("#define $mpidef");
ph("#endif");

for my $f (@ARGV) {
  proc_file($f);
}
ph("#endif");

open(G,">sct.h") || die "cannot open sct.h";
print G @h;
close(G);

