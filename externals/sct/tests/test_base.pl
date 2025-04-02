use strict;
use warnings;

my %ref_mpi_omp=();
my %ref_mpi=();
my %ref_serial=();
my @report_stack=();
my $ref_mpi_size = 0;
my $ref_omp_size = 0;
my $err_msg = "";
my $err = 0.0;
my $verbose = 1;
my $err_tol = 0.001;

sub run_prg {
  my @prg = @_;
  my @arg = split(' ', $_[0]);

  print "\n*** run_prg: prg = @prg ***\n";
  if ((-x @prg) || (-x $arg[0])) {
    my @out = `@prg`;
    if ($?) {
      my $exitcode = $? >> 8;
      die "ERROR: non-zero exit code ($exitcode) received from @prg"
    }
    print "run_prg: [out]=[@out]\n";
    return @out;
  } else {
    die "executable script expected in \@prg or \$arg[0]";
  }
}

sub parser_start {
  my %r;
  $r{thread_choice} = "";
  $r{nprocs} = 1;
  $r{nthreads} = 1;
  $r{proc_choice} = "";
  $r{sp_merging} = "";
  $r{cstats} = 0; #call-stats
  $r{proc_stats} = 0; #mpi parallel stats
  $r{thread_stats} = 0; #omp parallel stats
  $r{proc_details} = 0;
  $r{thread_details} = 0;
  $r{with_lbe} = 0;
  return %r;
}

sub parse {
  my @out = @_;

  %ref_mpi_omp=();
  %ref_mpi=();
  %ref_serial=();
  $ref_mpi_size = 0;
  $ref_omp_size = 0;
  @report_stack=();
  my %r=parser_start(); #report
  my $report_state = 0;
  my $stat_count = 0;
 parse_loop:for my $line (@out) {
    chomp($line);
    my $x =  $line . " ";
    print "parse:[$line]\n" if $verbose;

    if ($x =~ /^\s*(debug:|#)/) {
      next parse_loop;
    }

    if ($x =~ /^\s*ERROR/) {
        die "found error message";
    }

    if ($x =~ /^\s*ref_mpi_omp_dt:/) {
      if ($x =~ /^\s*ref_mpi_omp_dt:\s*label=\s*([^:]+?)\s*,\s*pid=\s*(\d+)\s*,\s*tid=\s*(\d+)\s*,\s*dt=\s*(\S+)/) {
        #print "mpi_omp: $1, $2, $3, $4\n";
        my ($label, $pid, $tid, $dt) = ($1, $2, $3, $4);
        push @{$ref_mpi_omp{$label}[$pid][$tid]},$dt;
        next parse_loop;
      } else {
        die "cannot parse $x\n";
      }
    }

    if ($x =~ /^\s*ref_mpi_dt:/) {
      if ($x =~ /^\s*ref_mpi_dt:\s*label=\s*([^:]+?)\s*,\s*pid=\s*(\d+)\s*,\s*dt=\s*(\S+)/) {
        my ($label, $pid, $dt) = ($1, $2, $3);
        push @{$ref_mpi{$label}[$pid]},$dt;
        #print "mpi: $1, $2, $3\n";
        next parse_loop;
      } else {
        die "cannot parse $x\n";
      }
    }

    if ($x =~ /^\s*ref_omp_dt:/) {
      if ($x =~ /^\s*ref_omp_dt:\s*label=\s*([^:]+?)\s*,\s*tid=\s*(\d+)\s*,\s*dt=\s*(\S+)/) {
        my ($label, $tid, $dt) = ($1, $2, $3);
        push @{$ref_mpi_omp{$label}[0][$tid]},$dt;
        print "omp: $1, $2, $3\n" if $verbose;
        next parse_loop;
      } else {
        die "cannot parse $x\n";
      }
    }

    if ($x =~ /^\s*ref_serial_dt:/) {
      if ($x =~ /^\s*ref_serial_dt:\s*label=\s*([^:]+?)\s*,\s*dt=\s*(\S+)/) {
        my ($label, $dt) = ($1, $2);
        push @{$ref_serial{$label}},$dt;
        print "serial: $1, $2\n" if $verbose;
        next parse_loop;
      } else {
        die "cannot parse $x\n";
      }
    }

    if ($x =~ /^\s*--> sct_report .*:/) {
      $report_state=1;
      print "(1) report_state=$report_state\n" if $verbose;
      $r{cstats} = 0;
      $r{proc_stats} = 0;
      $r{thread_stats} = 0;
      $r{proc_details} = 0;
      $r{thread_details} = 0;
      next parse_loop;
    }

    if ($report_state < 1) {
      print "ignore0:[$x]\n" if $verbose;
      next parse_loop;
    }

    if ($report_state == 1) {
      if ($x =~ /^\s*name\s*\|/) {
        $r{proc_details} = 1 if $x =~ /\sproc\s/;
        $r{thread_details} = 1 if $x =~ /\sthread\s/;
        $r{with_lbe} = 1 if $x =~ /\slbe/;
        $r{cstats} = 1 if $x =~ /\s\#calls\s/;
        $r{proc_details} = 1 if $r{proc_choice} eq "SCT_SELECT_ALL";
        $r{proc_stats} = 1 if $r{proc_choice} eq "SCT_REDUCE_ALL";
        $r{thread_details} = 1 if $r{thread_choice} eq "SCT_SELECT_ALL";
        $r{thread_stats} = 1 if $r{thread_choice} eq "SCT_REDUCE_ALL";
        $r{zero_stats} = 1 if not ($r{cstats} || $r{proc_stats} || $r{thread_stats});
        $r{with_lbe} = 1 if ($r{proc_stats} || $r{thread_stats});
        print "proc_details=$r{proc_details}, thread_details=$r{thread_details}, cstats=$r{cstats}, with_lbe=$r{with_lbe}\n" if $verbose;
        next parse_loop;
      }
    }

    if ($x =~ /^\s*---/) {
      $report_state = ($report_state+1) % 3;
      print "(2) report_state=$report_state\n" if $verbose;
      if (!$report_state) {
        my %cp_r = %r;
        push @report_stack, \%cp_r;
        %r=parser_start();
      }
      next parse_loop;
    }

    if ($x =~ /^\s*(\S+)\s*=\s*(\S+)\s*$/) {
      $r{$1}=$2;
      print "set: [$1] to [$2]\n" if $verbose;
      next parse_loop;
    }

    if ($report_state < 2) {
      print "ignore1:[$x]\n" if $verbose;
      next parse_loop;
    }

    if ($x =~ s/^\s*([^\|]+?)\s*\|\s*(.+)$/$2/) {
      my $name = $1;
      $name =~ s/^\s*[L ]+//;
      my ($pid, $tid, $cnum) = (0, 0, 0);
      my ($min, $avg, $max, $sum) = (0, 0, 0, 0);
      my $lbe = 100;
      if ($r{proc_details}) {
        if ($x =~ s/^\s*(\d+)\s+//) {
          $pid=$1;
        } else {
          die "cannot parse proc_details in [line]=[$line]\n";
        }
      }
      if ($r{thread_details}) {
        if ($x =~ s/^\s*(\d+)\s+//) {
          $tid=$1;
        } else {
          die "cannot parse thread_details in [line]=[$line]\n";
        }
      }
      if ($r{proc_details} || $r{thread_details}) {
        die "missing separator" unless  $x =~ s/^\s*\|//;
      }

      if ($r{cstats}) {
        if ($x =~ s/^\s*(\d+)\s+//) {
          $cnum=$1;
        } else {
          die "cannot parse cstats in [line]=[$line]\n";
        }
        $stat_count = $cnum;
      } else {
        $stat_count = 1;
        $stat_count *= $r{nprocs} if $r{proc_stats};
        $stat_count *= $r{nthreads} if $r{thread_stats};
      }

      if ($r{cstats} || $r{proc_stats} || $r{thread_stats}) {
        if ($x =~ s/^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*//) {
          ($min, $avg, $max, $sum) = ($1, $2, $3, $4);
        } else {
          die "cannot parse x=[$x]\n";
        }
      } else {
        die "internal error" unless $r{zero_stats};
        if ($x =~ s/^\s*(\S+)\s*//) {
          ($min, $avg, $max, $sum) = (undef, undef, undef, $1);
        } else {
          die "cannot parse x=[$x]\n";
        }
      }
      if ($r{with_lbe}) {
        if ($x =~ s/^\s*(\S+)\s*//) {
          $lbe = $1;
        } else {
          die "cannot parse x=[$x]\n";
        }
      }
      $x =~ s/\s*//;
      die "incomplete parsing of [line]=$line" if $x;

      #print "stat_count = $stat_count\n";
      #print "name=[$name], pid=$pid, tid=$tid, cnum=$cnum, min=$min, avg=$avg, max=$max, sum=$sum, lbe=$lbe, x=[$x]\n";
      if ($r{zero_stats}) {
        #print "zero_stats: name=[$name], pid=$pid, tid=$tid, sum=$sum, x=[$x]\n";
        $r{rec}{$name}[$pid][$tid] = $sum;
      } else {
        $r{rec}{$name}[$pid][$tid] = [$stat_count, $min, $avg, $max, $sum, $lbe];
      }
      next parse_loop;
    }

    die "cannot parse [line]=[$line]\n";

  } # end of parse_loop


  $ref_mpi_size = 0;
  for my $k (keys %ref_mpi) {
    if (!$ref_mpi_size) {
      $ref_mpi_size = scalar @{$ref_mpi{$k}};
    } else {
      die "bad ref_mpi_size" unless $ref_mpi_size == scalar @{$ref_mpi{$k}};
    }
  }
  #print "ref_mpi_size=$ref_mpi_size\n";
  $ref_omp_size = 0;
  for my $k (keys %ref_mpi_omp) {
    if (!$ref_mpi_size) {
      $ref_mpi_size = scalar @{$ref_mpi_omp{$k}};
    } else {
      die "bad ref_mpi_size" unless $ref_mpi_size == scalar @{$ref_mpi_omp{$k}};
    }
    for (my $pid = 0; $pid < $ref_mpi_size; $pid++) {
      if (!$ref_omp_size) {
        $ref_omp_size = scalar @{$ref_mpi_omp{$k}[$pid]};
      } else {
        die "bad ref_mpi_size" unless $ref_omp_size == scalar @{$ref_mpi_omp{$k}[$pid]};
      }
    }
  }
  #print "ref_omp_sizes=$ref_omp_size\n";
  return 1;
}

sub gen_stats {
    my @tdata = @_;
    return (0, 0, 0, 0, 0, 100) unless @tdata;

    my $min = $tdata[0];
    my $max = $min;
    my $sum = 0;
    #print("num tdata=", scalar @tdata,"\n");
    for my $dt (@tdata) {
      $sum += $dt;
      $min = $dt if $dt < $min;
      $max = $dt if $dt > $max;
    }
    my $count = scalar(@tdata);
    my $avg = $sum / $count;
    my $lbe;
    if ($max>0.0) {
      $lbe= 100.0 * $avg / $max;
    } else {
      $lbe = 100.0;
    }
    return ($count, $min, $avg, $max, $sum, $lbe);
}

sub gen_sum {
    my @tdata = @_;
    return 0 unless @tdata;

    my $sum = 0;
    for my $dt (@tdata) {
      $sum += $dt;
    }
    return $sum;
}

sub gen_ref_report {
  my $r_pt = shift || 0;
  die "HASH ref expected" unless ref($r_pt) eq "HASH";
  my %r=%{$r_pt}; # guiding data for constructing reference data
  my $rec_pt = $r{rec} || 0;
  my %rec= %{$rec_pt};
  my %ref_merged = ();
  my %s = (); #this reference report;

  # serial times:
  for my $k (keys %ref_serial) {
    @{$ref_merged{$k}[0][0]} = @{$ref_serial{$k}};
    #print "ref_serial to ref_merged: k=$k, ref_merged = @{$ref_merged{$k}[0][0]}\n";
  }

  # mpi times:
  if ($ref_mpi_size) {
    if ( $r{program_mode} eq "MPI" or
         $r{sp_merging} eq "SCT_SP_SERIAL_ONLY" or
         $r{sp_merging} eq "SCT_SP_MERGE_SIMPLE"
       ) {
      for my $k (keys %ref_mpi) {
        for (my $pid = 0; $pid < $ref_mpi_size; $pid++) {
          if ($ref_mpi{$k}[$pid]) {
            @{$ref_merged{$k}[$pid][0]} = @{$ref_mpi{$k}[$pid]};
            #print "sp_merging1: k=$k, pid=$pid: ref_merged[0] = @{$ref_merged{$k}[$pid][0]}\n";
          }
        }
      }
    }
  }

  #omp times:
  if ($ref_omp_size) {
    if ($r{sp_merging} eq "SCT_SP_MERGE_SIMPLE"  ) {
      for my $k (keys %ref_merged) {
        for (my $pid = 0; $pid < $ref_mpi_size; $pid++) {
          for (my $tid = 1; $tid < $ref_omp_size; $tid++) {
            #print "zero ref_merged{$k}[0][$tid]\n";
            @{$ref_merged{$k}[$pid][$tid]} = () unless $ref_merged{$k}[$pid][$tid];
          }
        }
      }
    }

    if ( $r{sp_merging} eq "SCT_SP_PARALLEL_ONLY" or
         $r{sp_merging} eq "SCT_SP_MERGE_SIMPLE"  ) {
      #print "sp_merging parallel part\n" if $verbose;
      for my $k (keys %ref_mpi_omp) {
        $ref_merged{$k} = [] unless $ref_merged{$k};
        for (my $pid = 0; $pid < $ref_mpi_size; $pid++) {
          for (my $tid = 0; $tid < $ref_omp_size; $tid++) {
            @{$ref_merged{$k}[$pid][$tid]} = () unless $ref_merged{$k}[$pid][$tid];
            if ($ref_mpi_omp{$k}[$pid][$tid]) {
              push @{$ref_merged{$k}[$pid][$tid]}, @{$ref_mpi_omp{$k}[$pid][$tid]};
              #print "sp_merging2: k=$k, pid=$pid, tid=$tid: ref_merged = @{$ref_merged{$k}[$pid][$tid]}\n";
            }
          }
        }
      }
    }
  }


  # call stats
  if ($r{cstats}) {
    for my $k (keys %ref_merged) {
      for (my $pid = 0; $pid < scalar(@{$ref_merged{$k}}); $pid++) {
        for (my $tid = 0; $tid < scalar(@{$ref_merged{$k}[$pid]}); $tid++) {
          #$ref_mpi_cnum{$k}[$pid] = scalar( @{$ref_mpi{$k}[$pid]} );
          #print "ref_mpi_cnum{$k}[$pid] = $ref_mpi_cnum{$k}[$pid]\n";
          my ($stat_count, $min, $avg, $max, $sum, $lbe) =
            gen_stats(@{$ref_merged{$k}[$pid][$tid]});
          $lbe = 100 unless $r{with_lbe};
          $s{rec}{$k}[$pid][$tid] = [$stat_count, $min, $avg, $max, $sum, $lbe];
        }
      }
    }
    return \%s;
  }

  # call-sum:
  my %csum=();
  for my $k (keys %ref_merged) {
    for (my $pid = 0; $pid < scalar(@{$ref_merged{$k}}); $pid++) {
      for (my $tid = 0; $tid < scalar(@{$ref_merged{$k}[$pid]}); $tid++) {
        $csum{$k}[$pid][$tid] = gen_sum(@{$ref_merged{$k}[$pid][$tid]});
      }
    }
  }

  for my $k (keys %ref_merged) {
    my $np = scalar(@{$csum{$k}});
    my $nt = scalar(@{$csum{$k}[0]});
    for (my $pid = 0; $pid < $np; $pid++) {
      die "non uniform thread size" unless $nt == scalar(@{$csum{$k}[$pid]});
    }
    if ($r{thread_stats}) {
      if ($r{proc_stats}) {

        my @a = ();
        for (my $pid = 0; $pid < $np; $pid++) {
          for (my $tid = 0; $tid < $nt; $tid++) {
            push @a, $csum{$k}[$pid][$tid];
          }
        }
        my ($stat_count, $min, $avg, $max, $sum, $lbe) = gen_stats(@a);
        $s{rec}{$k}[0][0] = [$stat_count, $min, $avg, $max, $sum, $lbe];

      } else {

        for (my $pid = 0; $pid < $np; $pid++) {
          my @a = ();
          for (my $tid = 0; $tid < $nt; $tid++) {
            push @a, $csum{$k}[$pid][$tid];
          }
          my ($stat_count, $min, $avg, $max, $sum, $lbe) = gen_stats(@a);
          $s{rec}{$k}[$pid][0] = [$stat_count, $min, $avg, $max, $sum, $lbe];
        }

      }
    } else {
      if ($r{proc_stats}) {

        my @a = ();
        for (my $pid = 0; $pid < $np; $pid++) {
          for (my $tid = 0; $tid < $nt; $tid++) {
            push @{$a[$tid]}, $csum{$k}[$pid][$tid];
          }
        }
        for (my $tid = 0; $tid < $nt; $tid++) {
          my ($stat_count, $min, $avg, $max, $sum, $lbe) = gen_stats(@{$a[$tid]});
          $s{rec}{$k}[0][$tid] = [$stat_count, $min, $avg, $max, $sum, $lbe];
        }

      } else {

        die "unexpected reduction case" unless $r{zero_stats};
        for (my $pid = 0; $pid < $np; $pid++) {
          for (my $tid = 0; $tid < $nt; $tid++) {
            $s{rec}{$k}[$pid][$tid] = $csum{$k}[$pid][$tid];
          }
        }

      }
    }
  }

  return \%s;

  # dead code:

  # thread-reduction:
  my %tred=();
  if ($r{thread_stats}) {
    for my $k (keys %csum) {
      for (my $pid = 0; $pid < scalar(@{$csum{$k}}); $pid++) {
        $tred{$k}[$pid][0] = gen_sum(@{$csum{$k}[$pid]});
      }
    }
  } else {
    %tred = %csum;
  }

  # proc-sum:
  my %pred=();
  if ($r{proc_stats}) {
    for my $k (keys %tred) {
      my @tmp = ();
      for (my $pid = 0; $pid < scalar(@{$tred{$k}}); $pid++) {
        for (my $tid = 0; $tid < scalar(@{$tred{$k}[$pid]}); $tid++) {
          $tmp[$tid][$pid] = $tred{$k}[$pid][$tid];
        }
      }
      for (my $tid = 0; $tid < scalar(@tmp); $tid++) {
        $pred{$k}[0][$tid] = gen_sum(@{$tmp[$tid]});
      }
    }
  } else {
    %pred = %tred;
  }

  # reduction stats
  for my $k (keys %pred) {
    for (my $pid = 0; $pid < scalar(@{$pred{$k}}); $pid++) {
      for (my $tid = 0; $tid < scalar(@{$pred{$k}[$pid]}); $tid++) {
        print "reduction stats: k=$k, pid=$pid, tid=$tid\n";
        my ($stat_count, $min, $avg, $max, $sum, $lbe) =
          gen_stats(@{$pred{$k}[$pid][$tid]});
        $lbe = 100 unless $r{with_lbe};
        $s{rec}{$k}[$pid][$tid] = [$stat_count, $min, $avg, $max, $sum, $lbe];
      }
    }
  }
  return \%s;

}

sub check_err {
    my ($label, $val, $ref_val) = @_;
    my $eps = 0.000001;
    my $rel_e = abs($val-$ref_val)/($ref_val+$eps);
    my $abs_e = abs($val-$ref_val);
    my $e = $rel_e;
    $e = $abs_e if $abs_e < $rel_e;
    if ($e>$err_tol) {
      $err_msg .= "check_error($label, $val, $ref_val): $e\n";
    } else {
      print "check_err($label, $val, $ref_val): err=$e : okay\n" if $verbose;
    }
    return $e;
}

sub check_single_report {
  my $r_pt = shift;
  die "undefined argument" unless defined $r_pt;
  die "expected HASH ref" unless ref($r_pt) eq "HASH";
  print "check_single_report: r_pt=$r_pt\n" if $verbose;
  my %r = %{$r_pt};

  my $nprocs = $r{nprocs} || 1;
  my $nthreads = $r{nthreads} || 1;
  print "nprocs=$nprocs, nthreads=$nthreads\n" if $verbose;
  my $rec_pt = $r{rec} || 0;
  die "record reference missing" unless ref($rec_pt) eq "HASH";
  my %rec = %{$rec_pt};

  my $ref_r_pt = gen_ref_report($r_pt) || 0;
  die "reference missing" unless ref($ref_r_pt) eq "HASH";
  my %ref_r = %{$ref_r_pt};
  my $ref_rec_pt = $ref_r{rec} || 0;
  die "reference record reference missing" unless ref($ref_rec_pt) eq "HASH";
  my %ref_rec = %{$ref_rec_pt};

  for my $k (keys %rec) {
    print "check key: $k\n" if $verbose;
    #print "stat_count okay ($stat_count == $ref_stat_count)\n"
    for (my $pid = 0; $pid < scalar(@{$rec{$k}}); $pid++) {
      for (my $tid = 0; $tid < scalar(@{$rec{$k}[$pid]}); $tid++) {
        $err_msg = "";
        $err = 0.0;
        if ($r{zero_stats}) {
          #print "zero_stats: check $rec{$k}[$pid][$tid] with $ref_rec{$k}[$pid][$tid]\n";
          $err += check_err("value", $rec{$k}[$pid][$tid], $ref_rec{$k}[$pid][$tid]);
          if ($err > $err_tol) {
            print "err_msg: $err_msg\n";
            print "k=$k, pid=$pid, tid=$tid\n";
            print "rec: $rec{$k}[$pid][$tid]\n";
            print "ref: $ref_rec{$k}[$pid][$tid]\n";
            die "serious error count: $err\n";
          }
        } else {
          my $a_pt = $rec{$k}[$pid][$tid] || 0;
          my $ref_a_pt = $ref_rec{$k}[$pid][$tid] || 0;
          die "reference a_pt expected" unless ref($a_pt) eq "ARRAY";
          die "empty reference ref_a_pt data" unless ref($ref_a_pt) eq "ARRAY";
          my ($stat_count, $min, $avg, $max, $sum, $lbe) =
            @{$rec{$k}[$pid][$tid]};
          my ($ref_stat_count, $ref_min, $ref_avg, $ref_max, $ref_sum, $ref_lbe) =
            @{$ref_rec{$k}[$pid][$tid]};
          #$s{rec}{$k}[$pid][$tid] = [$stat_count, $min, $avg, $max, $sum, $lbe];

          if ($r{cstats}) {
            if ($stat_count != $ref_stat_count) {
              $err_msg .=  "stat_count: ($stat_count != $ref_stat_count)\n";
              $err += 1.0;
            } else {
              print "stat_count okay ($stat_count == $ref_stat_count)\n" if $verbose;
            }
          }
          $err += check_err("min", $min, $ref_min);
          $err += check_err("avg", $avg, $ref_avg);
          $err += check_err("max", $max, $ref_max);
          $err += check_err("sum", $sum, $ref_sum);
          $err += check_err("lbe", $lbe, $ref_lbe);
          #die "array ref expected" unless ref($a_pt) eq "ARRAY";
          #my [$cnum, $min, $avg, $max, $sum, $lbe]
          if ($err > $err_tol) {
            print "err_msg: $err_msg\n";
            print "k=$k, pid=$pid, tid=$tid\n";
            print "rec: ($stat_count, $min, $avg, $max, $sum, $lbe)\n";
            print "ref: ($ref_stat_count, $ref_min, $ref_avg, $ref_max, $ref_sum, $ref_lbe)\n";
            die "serious error count: $err\n";
          }

        }
      }
    }
  }

}

sub check_reports {
  my $n = scalar(@report_stack);
  print "check_reports: stacksize = $n\n" ;
  for my $r_ref (@report_stack) {
    check_single_report($r_ref);
  }
}

1;
