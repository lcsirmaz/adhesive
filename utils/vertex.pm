#########################################################################
#########################################################################

package vertex;
use strict;

##########################################################################
# this package checks if a given vector is an internal point of the
# projection. It does by trying to solve the LP for the  missing variables
# It uses the external program "gspx" with the first argument "-x", which,
# in turn, uses the glpk LP solver
#
# tha main routine is is_vertex($arr), which aborts if $arr is not internal
# 

sub val{
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

# compute the RHS for $sh, add if not all zero coeffs in $vars
sub add_SH_inequality {
   my($info,$sh,$arr)=@_;
   my $v=0; my $add=0;
   for my $i(1 .. 63){
      next if($sh->[$i]==0);
      if($info->{vars}->[$i]>=0){
         $add=1;
         die "is_vertex: var $i cannot be set\n" if($arr->[$i]!=0);
      } else { $v -= $sh->[$i]*$arr->[$i]; }
   }
   if($add){
       push @{$info->{M}},$sh; push @{$info->{RHS}},$v;
   } elsif($v>0){
       die "is_vertex: supplied point does not satisfy Shannon\n";
   }
}

# check if $arr->[1 .. 63] is an internal point
sub is_vertex {
   my($arr)=@_;
   # $vars->[$i]=-1 if fixed value, otherwise it grows from 0
   my $vars=getarr(); my $cols=0;
   for my $i(0 .. 63){
      $vars->[$i]=-1;
      if(($i&val("xy"))!=0 && ($i&val("uv"))!=0){
        $vars->[$i]=$cols; $cols++;
      }
   }
   #   die "is_vertex: $cols!=36\n"; ## debugging
   # create Shannon inequalities plus (xy,uv|ab)=0
   # M: Shannon ineuqlities, RHS: right hand side after >=0
   my $info={ vars=>$vars, M=>[], RHS=>[] };
   # inequalities (i,j|A)>=0
   for my $A(0 .. 63){
     for (my $i=1;$i<63;$i<<=1){
        next if($A&$i);
        for (my $j=$i; $j<63;$j<<=1){
           next if($A&$j); next if($i==$j);
           my $c=getarr(); $c->[$i|$A]=1; $c->[$j|$A]=1;
           $c->[$i|$j|$A]=-1; if($A){$c->[$A]=-1; }
           add_SH_inequality($info,$c,$arr);
        }
     }
   }
   # inequalities (i|M-i)>=0
   for(my $i=1;$i<63;$i<<=1){
      my $c=getarr(); $c->[63]=1; $c->[63^$i]=-1;
      add_SH_inequality($info,$c,$arr);
   }
   # (xy,uv|ab)=0, last row
   my $c=getarr(); $c->[val("xyab")]=1; $c->[val("uvab")]=1;
   $c->[val("xyuvab")]=-1; $c->[val("ab")]=-1;
   add_SH_inequality($info,$c,$arr);
   my $rows=scalar @{$info->{M}};
   my $fname=`mktemp -t -q vertexXXXXXX.spx`;
   chomp $fname;
   my $resname=$fname; $resname =~ s/\.spx/\.res/;
   open(SPX,">",$fname) || die "is_vertex: cannot open $fname\n";
     # cols, rows, V
     print SPX "$cols\n$rows\n0\n";
     # matrix M, first row, second row, etc
     for my $i(0 .. $rows-1){
       for my $j(1 .. 63){
         next if($info->{vars}->[$j]<0);
         print SPX $info->{M}->[$i]->[$j],"\n";
       }
     }
     # right hand side, >= except for the last row, where we have =
     for my $i(0 .. $rows-1){
        print SPX ($i==$rows-1?"":"g"),$info->{RHS}->[$i],"\n";
     }
   close(SPX);
   system("gspx","-x",$fname,$resname);\
   open(RES,"<",$resname)||die "is_vertex: no result for LP $fname\n";
     my $yes=(<RES> =~ /V=XXXX/) ? 0 : 1;
   close(RES);
   unlink $fname,$resname;
   return $yes;
}

# indicate that all is fine

1;

__END__
