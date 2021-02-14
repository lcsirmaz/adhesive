#!/usr/bin/perl -W -I ../utils
##
## (xy,uv|ab)=0
##  add: x|aby, y|abx, u|abv, v|abu
##  
##  (xy,uv|ab)=0 decreases the number of variables
##

use strict;

my $case=1;      # 1 or 2, giving a|xyb or a|uvb; b|axy or b|auv
my $ingleton1=1; # 1,.., 6 which Ingleton on abxy is <=0
my $ingleton2=1; # 1,.., 6 which Ingleton on abuv is <=0

my $facetEqFile = "../allfacets.txt";
my $polcoJar = "../utils/polco.jar";

sub val {
   my $txt=shift;
   my $v=0;
   $txt =~ s/\s//g;
   while($txt){
     if($txt !~ /^[abxyuv]/){ die "val: wrong argument $txt\n"; }
     $v|=1 if ($txt =~ s/^a//);
     $v|=2 if ($txt =~ s/^b//);
     $v|=4 if ($txt =~ s/^x//);
     $v|=8 if ($txt =~ s/^y//);
     $v|=16 if($txt =~ s/^u//);
     $v|=32 if($txt =~ s/^v//);
   }
   return $v;
}

sub name{
   my $v=shift;
   my $txt="";
   $txt .= "a" if($v&1);
   $txt .= "b" if($v&2);
   $txt .= "x" if($v&4);
   $txt .= "y" if($v&8);
   $txt .= "u" if($v&16);
   $txt .= "v" if($v&32);
   return $txt;
}

sub getarr { # return an empty array of size 64
   my @c=(0)x64;
   return \@c;
}

sub mkhash { # create a string from $arr
   my ($arr)=@_; 
   my $h="";
   foreach my $v(@$arr){
      $h.=",";
      next if $v==0;
      $h.=$v."";
   }
   return $h;
}

##########################################################
# We maintain a replacement table as follows:
#   $info->{T}->{val("old")} = { val("new") => 1 }

sub translate { # replace $c[0..63] with the defined variables
    my($info,$c)=@_;
    my($T,$cc)=($info->{T},$info->{cc});
    if(! defined $cc){$cc=getarr(); $info->{cc}=$cc; }
    for my $i(0 .. 63){$cc->[$i]=$c->[$i]; $c->[$i]=0;}
    for my $k(0 .. 63){
       my $v=$cc->[$k];
       next if($v==0);
       if(! defined $T->{$k}){ $c->[$k]+=$v; }
       else{ foreach my $t(keys %{$T->{$k}}){
           $c->[$t]+= $v*$T->{$k}->{$t};
       }}
    }
    return $c;
}

###################################################################
# first stage: add x|aby, etc. It gives equal values for several
#  variables.
# define $a and $b to be equal where $a is a superset of $b
#  it implies ax = bx as 0<=bx-ax<=b-a=0
sub make_eqvalue {
   my($info,$a,$b)=@_;
   $a=val($a); $b=val($b);
   die "make_eqvalue: $a is not a superset of $b\n" if(($a&$b)!=$b);
   for my $i(1 .. 63){
      my($aa,$bb)=($a|$i,$b|$i);
      next if($aa==$bb);
      if(defined $info->{T}->{$bb}){
         $bb = (keys %{$info->{T}->{$bb}})[0];
         next if($aa==$bb);
      }
      if(defined $info->{T}->{$aa}){
         $aa = (keys %{$info->{T}->{$aa}})[0];
         next if($aa==$bb);
      }
      # replace $aa by $bb in the definitions
      foreach my $t(keys %{$info->{T}}){
         if($aa == (keys %{$info->{T}->{$t}})[0]){
             $info->{T}->{$t}={$bb => 1};
         }
      }
      $info->{T}->{$aa} = { $bb => 1 };
   }
}

# generating Shannon inequalities for (xyab) and (abuv)
#   $mask is either val("xyab") or val("abuv")
sub genSH {
    my($info,$mask)=@_;
    ## inequalities (i,j|A)>=0
    for my $A(0 .. 63){
       next if(($A&$mask)!=$A);
       for(my $i=1;$i<63;$i<<=1){
          next if(($A&$i)!=0);
          for(my $j=$i; $j<63; $j<<=1){
          next if($i==$j);
          next if(($A&$j)!=0);
          my $Aij=$A|$i|$j; next if(($Aij&$mask)!=$Aij);
          my $c=getarr();
          $c->[$A]=-1; $c->[$A|$i]=1; $c->[$A|$j]=1; $c->[$Aij]=-1;
          $c->[0]=0;
          translate($info,$c);
          my $h=mkhash($c);
          next if($info->{hash}->{$h});
          $info->{hash}->{$h}=1;
          push @{$info->{A}},$c;
       }}
    }
    ## inequalities (i|M-i)>=0
    for(my $i=1;$i<63;$i<<=1){
      next if(($i&$mask)==0);
      my $c=getarr();
      $c->[$mask]=1; $c->[$mask^$i]=-1;
      translate($info,$c);
      my $h=mkhash($c);
      next if($info->{hash}->{$h});
      $info->{hash}->{$h}=1;
      push @{$info->{A}},$c;
    }    
}

## generating matrix A: inequalities with zero coeffs
sub genA {
    my($info)=@_;
    $info->{hash}->{mkhash(getarr())}=1; # all zero
    make_eqvalue($info,"abxy","abx");
    make_eqvalue($info,"abxy","aby");
    make_eqvalue($info,"abuv","abu");
    make_eqvalue($info,"abuv","abv");
    if($case==1){ ##  a|xyb, b|axy
       make_eqvalue($info,"abxy","bxy");
       make_eqvalue($info,"abxy","axy");
    } elsif($case==2){ ## a|xyb, b|auv
       make_eqvalue($info,"abxy","bxy");
       make_eqvalue($info,"abuv","auv");
    } else {
       die "genA: wrong value for case=$case\n";
    }
    genSH($info,val("abxy"));
    genSH($info,val("abuv"));
}

##############################################################################
# add Ingleton inequality

sub add3 {  ## add n*(a,b|c) fo $cc
    my($cc,$a,$b,$c,$n)=@_;
    $cc->[$a|$b|$c] -= $n;
    $cc->[$a|$c] += $n;
    $cc->[$b|$c] += $n;
    if($c){$cc->[$c]-= $n;}
}
# add n*[a,b,c,d]>=0 to {A}
sub add_Ingleton {
    my($info,$a,$b,$c,$d,$n)=@_;
    my $cc=getarr();
    $a=val($a); $b=val($b); $c=val($c); $d=val($d);
    add3($cc,$a,$b,0,-$n); # -(a,b)+(a,b|c)+(a,b|d)+(c,d)
    add3($cc,$a,$b,$c,$n);
    add3($cc,$a,$b,$d,$n);
    add3($cc,$c,$d,0,$n);
    translate($info,$cc);
    push @{$info->{A}},$cc;
}

#############################################################################
# read known facets
#  the are stored as this:
#    -x-u+uv-a+2xa+ua-va-2b+2xb+ub+vb-uvb+3ab-3xab-uab
sub txttoarr {  # string to an array
   my($txt)=@_;
   $txt =~ s/\s//g;
   my @arr=(0)x64;
   my $n=1; my $sign=1;
   while($txt){
      if($txt =~ s/^\-//){ $sign=-1; next; }
      if($txt =~ s/^\+//){ $sign=1; next; }
      if($txt =~ s/^(\d+)//){ $n=$1; next; }
      if($txt =~ s/^([abuvxy]+)//){ my $idx=val($1);
        $arr[$idx]+=$sign*$n; $sign=1; $n=1;
        next; }
      die "syntax error in expression $txt\n";
   }
   return \@arr;
}

sub read_faceteqs {
    my($info,$file)=@_;
    open(FILE,"<",$file) || die "Cannot open facet file $file\n";
     while(<FILE>){
       chomp;
       my $c=txttoarr($_); translate($info,$c);
       my $h=mkhash($c);
       next if(defined $info->{hash}->{$h});
       $info->{hash}->{$h}=1;
       push @{$info->{A}},$c;
     }
    close(FILE);
}

################################################################################
#  create var indices for matrix A
#  {inA}->[$i]=-1 if not used, otherwise the serial number
sub create_varindex {
    my($info)=@_;
    $info->{inA}=[];
    for my $i(0 .. 63){$info->{inA}->[$i]=-1; }
    foreach my $eq(@{$info->{A}}){
      foreach my $i(1 .. 63){
         next if($eq->[$i]==0);
         $info->{inA}->[$i]=0;
      }
    }
    my $cnt=0;
    for my $i(1 .. 63){
       next if($info->{inA}->[$i]<0);
       $info->{inA}->[$i]=$cnt; $cnt++;
    }
    $info->{ninA}=$cnt;
}

##################################################################################
# filter out redundant inequalities
#  if $info->{i}->[$j]==1 it is a consequence of the others
sub try_indep {
    my($info,$idx)=@_;
    my $rows=$info->{ninA};
    my $cols=scalar @{$info->{A}};
    if(!defined $info->{i}){ $info->{i}=[]; for my $i(0 .. $cols-1){$info->{i}->[$i]=0;}}
    my $cnt=0; for my $i(0 .. $cols-1){ $cnt++ if($info->{i}->[$i]==0); }
    my $fname=`mktemp -t -q indXXXXXX.spx`;
    chomp $fname;
    my $resname=$fname; $resname =~ s/\.spx$/\.res/;
    open(SPX,">",$fname) || die "Cannot open $fname for writing\n";
      # cols, rows, V
      print SPX "",$cnt-1,"\n",$rows,"\n",0,"\n";
      # Matrix first row, second row, ... last row
      for my $j(0 .. 63){
         next if($info->{inA}->[$j]<0);
         for my $col(0 .. $cols-1){
            next if($col==$idx); ## this is the RHS
            next if($info->{i}->[$col]!=0);
            my $eq=$info->{A}->[$col];
            print SPX $eq->[$j]>=0?"  ":" ",$eq->[$j];
         }
         print SPX "\n";
      }
      # right hand side:  =
      for my $j(0 .. 63){
         next if($info->{inA}->[$j]<0);
         my $eq=$info->{A}->[$idx];
         print SPX $eq->[$j],"\n";
      }
    close(SPX);
    system("gspx","-x",$fname,$resname);
    open(RES,"<",$resname)|| die "no result for LP $fname\n";
      my $yes=(<RES> =~ /V=XXXX/) ? 0 : 1;
    close(RES);
    unlink $fname,$resname;
    if($yes){$info->{i}->[$idx]=1; }
}
##################################################################################
# create all inequalities for <case><ingleton> as a two-digit number

sub create_inequalities {
    my($arg)=@_;
    if($arg !~ /^([1-2])([1-6])([1-6])$/){ die "create_inequalities: wrong arg $arg\n"; }
    $case=$1; $ingleton1=$2; $ingleton2=$3;
    my $info={ A=>[], hash=>{}, T=>{} };
    genA($info);
    ## add Ingleton inequalities for xyab
    if($ingleton1==1){ add_Ingleton($info,"a","b","x","y",-1);}
    elsif($ingleton1==2){ add_Ingleton($info,"a","x","b","y",-1);}
    elsif($ingleton1==3){ add_Ingleton($info,"a","y","b","x",-1);}
    elsif($ingleton1==4){ add_Ingleton($info,"b","x","a","y",-1);}
    elsif($ingleton1==5){ add_Ingleton($info,"b","y","a","x",-1);}
    elsif($ingleton1==6){ add_Ingleton($info,"x","y","a","b",-1);}
    else {die "ingleton=$ingleton1 is not 1..6\n"; }
    ## add Ingleton inequalities for abuv
    if($ingleton2==1){ add_Ingleton($info,"a","b","u","v",-1);}
    elsif($ingleton2==2){ add_Ingleton($info,"a","u","b","v",-1);}
    elsif($ingleton2==3){ add_Ingleton($info,"a","v","b","u",-1);}
    elsif($ingleton2==4){ add_Ingleton($info,"b","u","a","v",-1);}
    elsif($ingleton2==5){ add_Ingleton($info,"b","v","a","u",-1);}
    elsif($ingleton2==6){ add_Ingleton($info,"u","v","a","b",-1);}
    else {die "ingleton=$ingleton2 is not 1..6\n"; }
    # add facet equations
    read_faceteqs($info,$facetEqFile);
    $info->{hash}={}; # no more used
    create_varindex($info);
    # filter out dependent inequalities
    for my $i(0 .. -1+scalar @{$info->{A}}){ try_indep($info,$i); }
    return $info;
}

# swap last 60 lines to the beginning
sub print_MF {
    my($info,$file)=@_;
    open(FILE,">",$file) || die "print_MF: cannot create file $file\n";
    my $total=-1+scalar @{$info->{A}};
    die "print_MF: total=$total is <60\n" if($total<62);
    for my $i($total-60 .. $total){
        next if($info->{i}->[$i]!=0);
        for my $j(1 .. 63){
            next if($info->{inA}->[$j]<0);
            print FILE " ",$info->{A}->[$i]->[$j];
        }
        print FILE "\n";
    }
    for my $i(0 .. $total-61){
        next if($info->{i}->[$i]!=0);
        for my $j(1 .. 63){
            next if($info->{inA}->[$j]<0);
            print FILE " ",$info->{A}->[$i]->[$j];
        }
        print FILE "\n";
    }
    close(FILE);
}

################################################################################
# parse the result of POLCO
# it is a tab-separated list of  integers and fractions (all numbers are >=0)
# it is turned into an array of integers

sub gcd { # recursive procedure to compute the gcd of two >=0 integers
    my($u,$v)=@_;
    if($u==0 || $u==$v){ return $v; }
    if($v==0){ return $u; }
    if($v==1|| $u==1){ return 1; }
    if($u<$v){ return gcd( $v % $u, $u); }
    return gcd($u % $v, $v); # $v<$u
}

# from a tab-separated list create an array
#    3	8/3	13/3	7/3	4
sub parse_polco_line {
    my($line)=@_;
    my @a=split(/\s+/,$line); 
    my $den=1; # denominator
    foreach my $v(@a){
      next if($v !~ m#/(\d+)#);
      my $i=$1; my $j=gcd($den,$i);
      if($j==1){ $den *= $i; }
      elsif($j==$i){;}
      else{ $den = int(0.1+($den*$i)/($j+0.0)); }
    }
    foreach my $i(0 .. -1+scalar @a){
       if($a[$i] =~ m#(\d+)/(\d+)#){
         my($x,$y)=($1,$2);
         $a[$i]= $x* int(0.1+$den/$y);
       } else { $a[$i] *= $den; }
    }
    return \@a;
}

##################################################################################
# check if vertices returned by POLCO are indeed internal points
#
use vertex;

sub check_vertices {
    my($info,$file)=@_;
    my $c=getarr();
    open(PRES,"<",$file)|| die "Cannot open POLCO result $file\n";
      while(<PRES>){
        chomp;
        my $a=parse_polco_line($_);
        if(scalar @$a != $info->{ninA}){
           die "POLCO result line in $file has wrong dimension\n";
        }
        for my $j(0 .. 63){
            if($info->{inA}->[$j]<0) {$c->[$j]=0;}
            else{ $c->[$j]=$a->[$info->{inA}->[$j]]; }
        }
        # figure out the value of replaced values
        foreach my $k(keys %{$info->{T}}){
           foreach my $v(keys %{$info->{T}->{$k}}){
               $c->[$k] += $info->{T}->{$k}->{$v}*$c->[$v]
           }
        }
        # finally call is_vertex()
        my $yes=vertex::is_vertex($c);
        die "Not a vertex\n" if($yes==0);
      }
    close(PRES);
}

#################################################################################
# and the main routine
#

sub run_all {
    for my $x(1 .. 2){for my $y(1 .. 6){for my $z(1 .. 6){
       my $info=create_inequalities("$x$y$z");
       my $file="/tmp/FF$x$y$z";
       print_MF($info,"$file.dat");
       system("java","-Xms3g","-jar",$polcoJar,"-level","OFF",
           "-kind","text","-iq","$file.dat","-out","text","$file.res");
       check_vertices($info,"$file.res");
       unlink "$file.dat","$file.res";
    }}}
    print "*** DONE ***\n";
    exit 0;
}

if(scalar @ARGV>0){ $facetEqFile=$ARGV[0];}
run_all();

__END__

