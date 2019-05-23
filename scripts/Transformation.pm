#!/usr//bin/perl -w

package Transformation;

use strict;

sub multiply_transforms {
  my ($q1_0, $q1_1, $q1_2, $q1_3) = xyz_to_quaternion($_[0], $_[1], $_[2]);
  my ($q2_0, $q2_1, $q2_2, $q2_3) = xyz_to_quaternion($_[6], $_[7], $_[8]);

  my $q_0 = $q1_0*$q2_0 - $q1_1*$q2_1 - $q1_2*$q2_2 - $q1_3*$q2_3;
  my $q_1 = $q1_0*$q2_1 + $q1_1*$q2_0 + $q1_2*$q2_3 - $q1_3*$q2_2;
  my $q_2 = $q1_0*$q2_2 - $q1_1*$q2_3 + $q1_2*$q2_0 + $q1_3*$q2_1;
  my $q_3 = $q1_0*$q2_3 + $q1_1*$q2_2 - $q1_2*$q2_1 + $q1_3*$q2_0;

  my ($rx, $ry, $rz) = rotation_angles($q_0, $q_1, $q_2, $q_3);
  my ($t0, $t1, $t2) = applyRotation($q1_0, $q1_1, $q1_2, $q1_3, $_[9], $_[10], $_[11]);
  $t0 += $_[3];
  $t1 += $_[4];
  $t2 += $_[5];

  return ($rx, $ry, $rz, $t0, $t1, $t2);
}

sub xyz_to_quaternion {
  my $xr = $_[0]; my $yr = $_[1]; my $zr = $_[2];

  my $cx = cos($xr); my $cy = cos($yr); my $cz = cos($zr);
  my $sx = sin($xr); my $sy = sin($yr); my $sz = sin($zr);
  my $m00 = $cz*$cy;
  my $m11 = -$sy*$sx*$sz + $cx*$cz;
  my $m22 = $cy*$cx;

  my $q0 = sqrt(1+$m00+$m11+$m22)/2.0;
  my $q1 = sqrt(1+$m00-$m11-$m22)/2.0;
  my $q2 = sqrt(1-$m00+$m11-$m22)/2.0;
  my $q3 = sqrt(1-$m00-$m11+$m22)/2.0;

  if($cy*$sx + $sy*$cx*$sz + $sx*$cz < 0.0) { $q1 = -$q1; }
  if($sz*$sx - $sy*$cx*$cz - $sy < 0.0)     { $q2 = -$q2; }
  if($sz*$cy + $sy*$sx*$cz + $sz*$cx < 0.0) { $q3 = -$q3; }

  return ($q0, $q1, $q2, $q3);
}

sub rotation_angles {
  my $q0 = $_[0]; my $q1 = $_[1]; my $q2 = $_[2]; my $q3 = $_[3];

  my $matrix32 = 2*($q2*$q3 + $q0*$q1);
  my $matrix33 = $q0*$q0 - $q1*$q1 - $q2*$q2 + $q3*$q3;
  my $matrix31 = 2*($q1*$q3 - $q0*$q2);
  my $matrix21 = 2*($q1*$q2 + $q0*$q3);
  my $matrix11 = $q0*$q0 + $q1*$q1 - $q2*$q2 - $q3*$q3;

  my $x = atan2($matrix32, $matrix33);
  my $y = atan2($matrix31, sqrt($matrix21*$matrix21 + $matrix11*$matrix11));
  my $z = atan2($matrix21, $matrix11);
  return ($x, $y, $z);
}

sub applyRotation {
  my $q0 = $_[0]; my $q1 = $_[1]; my $q2 = $_[2]; my $q3 = $_[3];
  my $v0 = $_[4]; my $v1 = $_[5]; my $v2 = $_[6];

  my $tv0 = ($q0*$q0+$q1*$q1-$q2*$q2-$q3*$q3)*$v0 + 2*($q1*$q2-$q0*$q3)*$v1 + 2*($q1*$q3+$q0*$q2)*$v2;
  my $tv1 = 2*($q1*$q2+$q0*$q3)*$v0 + ($q0*$q0-$q1*$q1+$q2*$q2-$q3*$q3)*$v1 + 2*($q2*$q3-$q0*$q1)*$v2;
  my $tv2 = 2*($q1*$q3-$q0*$q2)*$v0 + 2*($q2*$q3+$q0*$q1)*$v1 + ($q0*$q0-$q1*$q1-$q2*$q2+$q3*$q3)*$v2;

  return ($tv0, $tv1, $tv2);
}

1;
