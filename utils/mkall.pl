#!/usr/bin/perl -w
use strict;

## this program reads "conditions.txt" and generates all symmteric versions
## a<->b, x<->y, u<->v, uv<->xy
##

##############################################################################
##############################################################################
## each (non-empty) subset of xyuvab is a variable; they are numbered 0 .. 63
##
## vv("text") -> numerical; vname(num) -> text

sub vval { 
  my($txt)=@_; my $v=0;
    $v|=1 if($txt =~ /x/);
    $v|=2 if($txt =~ /y/);
    $v|=4 if($txt =~ /u/);
    $v|=8 if($txt =~ /v/);
    $v|=16 if($txt =~ /a/);
    $v|=32 if($txt =~ /b/);
    return $v;
}
sub vname {
   my($v)=@_; my $txt="";
   $txt .="x" if($v&1);
   $txt .="y" if($v&2);
   $txt .="u" if($v&4);
   $txt .="v" if($v&8);
   $txt .="a" if($v&16);
   $txt .="b" if($v&32);
   return $txt;
}

# in string $txt swap chars $x and $y
sub swaptxt {
    my($x,$y,$txt)=@_;
    $txt =~ s/$x/A/g;
    $txt =~ s/$y/B/g;
    $txt =~ s/A/$y/ge;
    $txt =~ s/B/$x/ge;
    return vname(vval($txt));
}
# swap $x and $y in a text, return a translation array
sub swaparr { # $x<->$y
    my($x,$y,$arr)=@_; my $s=[];
    for my $i(0..63){
        $s->[vval(swaptxt($x,$y,vname($i)))]=$arr->[$i];
    }
    return $s;
}
sub swappair { # xy <-> uv
    my($arr)=@_; my $s=[];
    for my $i(0 .. 63){
       my $sw=swaptxt("y","v",swaptxt("x","u",vname($i)));
       $s->[vval($sw)]=$arr->[$i];
    }
    return $s;
}
##################################################################################
##################################################################################
# converting an expression to an array, and an array to an expression
sub arrtotxt {
   my($arr)=@_;
   my $txt="";
   for my $i(1 .. 63){
     my $v=$arr->[$i];
     next if($v==0);
     if($v==1){$txt .="+";}
     elsif($v==-1){$txt.="-";}
     elsif($v>0){$txt .="+$v";}
     else{$txt .= "$v"; }
     $txt .= vname($i);
   }
   $txt =~ s/^\+//;
   return $txt;
}
# add $n*($A,$B|$C) to $arr[]
sub add3toarr {
    my($arr,$A,$B,$C,$n)=@_;
    $arr->[vval("$A$C")]+=$n; $arr->[vval("$B$C")]+=$n;
    $arr->[vval("$A$B$C")]-=$n;
    if($C){$arr->[vval("$C")]-=$n;}
}
sub exprtoarr {
    my ($txt)=@_; $txt =~ s/\s//g;
    my @arr=(0)x64;
    my $n=1; 
    while($txt){
      $n=1;
      $txt =~ s/^\+//;
      if($txt =~ s/^(\d+)//){ $n=$1; }
      if($txt =~ s/^\[([abuvxy]),([abuvxy]),([abuvxy]),([abuvxy])\]//){
         my ($a,$b,$c,$d)=($1,$2,$3,$4);
         add3toarr(\@arr,$a,$b,"",-$n);
         add3toarr(\@arr,$a,$b,$c,$n);
         add3toarr(\@arr,$a,$b,$d,$n);
         add3toarr(\@arr,$c,$d,"",$n);
         next;
      }
      if($txt =~ s/^\(([abuvxy]+),([abuvxy]+)\|([abuvxy]+)\)//){
         my($a,$b,$c)=($1,$2,$3);
         add3toarr(\@arr,$a,$b,$c,$n);
         next;
      }
      die "wrong expression: $txt\n";
    }
    return \@arr;
}
#####################################################################################
#####################################################################################
# add an inequality and all symmetric version; check if one gets a new one
#  $C->{text}=serial_number
#
sub add_inequality {
   my($C,$n,$txt)=@_;
   my $a0=exprtoarr($txt);
    my $p0=arrtotxt($a0);
    if(defined $C->{$p0}){ die "ineq $txt is repeated as ($C->{$p0})\n"; }
    $C->{$p0}=$n;
    my $a1=swaparr("x","y",$a0);
    my $a2=swaparr("u","v",$a0);
    my $a3=swaparr("u","v",$a1);
    my $a4=swaparr("a","b",$a0);
    my $a5=swaparr("a","b",$a1);
    my $a6=swaparr("a","b",$a2);
    my $a7=swaparr("a","b",$a3);
    foreach my $a($a1,$a2,$a3,$a4,$a5,$a6,$a7){
      my $p=arrtotxt($a);
      if(defined $C->{$p} && $C->{$p}!=$n){
         die "symmetrical ineq $p again\n";
      }
      $C->{$p}=$n;
    }
    $a0=swappair($a0);
    $a1=swaparr("x","y",$a0);
    $a2=swaparr("u","v",$a0);
    $a3=swaparr("u","v",$a1);
    $a4=swaparr("a","b",$a0);
    $a5=swaparr("a","b",$a1);
    $a6=swaparr("a","b",$a2);
    $a7=swaparr("a","b",$a3);
    foreach my $a($a0,$a1,$a2,$a3,$a4,$a5,$a6,$a7){
      my $p=arrtotxt($a);
      if(defined $C->{$p} && $C->{$p}!=$n){
         die "symmetrical ineq $p again\n";
      }
      $C->{$p}=$n;
    }
}
#########################################################################################
#########################################################################################
# read inequalities, and then print out the full list
sub read_ineq_file {
    my($fname)=@_;
    my $C={}; my $n=1;
    open(FILE,"<".$fname)||die "Cannot open $fname\n";
    while(<FILE>){
       chomp;
       next if(/^$/);
       next if(/^#/);
       add_inequality($C,$n,$_); $n++;
    }
    close(FILE);
    return $C;
}
sub report_all {
    my($C)=@_;
    foreach my $k(sort {length($a)<=>length($b) || $a cmp $b} keys %$C){
        print "$k\n";
    }
}

#############################################################################################
#############################################################################################
# and execute ...

my $C=read_ineq_file("conditions.txt"); report_all($C); exit 3;



__END__




