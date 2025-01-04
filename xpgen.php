<?php

$MOD = gmp_init("0x16A6B036D7F2A79");

// plus 1 for comparassion
$max_uint64 = gmp_init("0xFFFFFFFFFFFFFFFF");
$max_uint64 = gmp_add($max_uint64,1); 

$max_uint32 = gmp_init("0xFFFFFFFF");
$max_uint32 = gmp_add($max_uint32,1); 

$max_uint16 = gmp_init("0xFFFF");
$max_uint16 = gmp_add($max_uint16,1); 

$max_uint8 = gmp_init("0xFF");
$max_uint8 = gmp_add($max_uint8,1); 

function gmp_shiftl($x,$n) { // shift left
    if(gettype($n)=="object") { // its a gmp number
        $n = gmp_intval($n);
    }
    return(gmp_mul($x,gmp_pow(2,$n)));
}
  
function gmp_shiftr($x,$n) { // shift right
    if(gettype($n)=="object") {
        $n = gmp_intval($n);
    }
    return(gmp_div($x,gmp_pow(2,$n)));
}


function uint64($in) { // makes number overflow like uint64
    global $max_uint64;
    
    $out = gmp_mod($in, $max_uint64);

    // underflow
    if (gmp_cmp($out, gmp_init(0)) < 0) {
        $out = gmp_add($out, $max_uint64);
    }


    return $out;
}

function uint32($in) { // makes number overflow like uint64
    global $max_uint32;
    
    $out = gmp_mod($in, $max_uint32);

    // underflow
    if (gmp_cmp($out, gmp_init(0)) < 0) {
        $out = gmp_add($out, $max_uint32);
    }


    return $out;
}

$max_int64 = gmp_init("9223372036854775807");
$max_int64 = gmp_add($max_int64,1); 
$min_int64 = gmp_init("-9223372036854775808");

function int64($in) {
    global $max_int64, $min_int64;
    
    // Handle overflow
    if (gmp_cmp($in, $max_int64) >= 0) {
        return gmp_sub($in, gmp_mul($max_int64, gmp_init(2)));
    }
    // Handle underflow
    elseif (gmp_cmp($in, $min_int64) < 0) {
        return gmp_add($in, gmp_mul($max_int64, gmp_init(2)));
    }
    return $in;
}

function uint8($in) {
    global $max_uint8;
    
    $out = gmp_mod($in, $max_uint8);

    // underflow
    if (gmp_cmp($out, gmp_init(0)) < 0) {
        $out = gmp_add($out, $max_uint8);
    }


    return $out;
}

function uint16($in) {
    global $max_uint16;
    
    $out = gmp_mod($in, $max_uint16);

    // underflow
    if (gmp_cmp($out, gmp_init(0)) < 0) {
        $out = gmp_add($out, $max_uint16);
    }


    return $out;
}


function residue_add($x, $y) {
    global $MOD;

	$z = uint64(gmp_add($x, $y));

    if(gmp_cmp($z, $MOD)>=0) {
        $z = uint64(gmp_sub($z, $MOD));
    }

	return $z;
}

function residue_sub($x, $y) {
    global $MOD;

	$z = uint64(gmp_sub($x, $y));

    if(gmp_cmp($x, $y)<0) {
        $z = uint64(gmp_add($z, $MOD));
    }
    
	return $z;
}

function __emulu($x, $y) {
    return uint64(gmp_mul($x,$y));
}

function __umul128($multiplier, $multiplicand, &$product_hi)  {
	// multiplier   = ab = a * 2^32 + b
	// multiplicand = cd = c * 2^32 + d
	// ab * cd = a * c * 2^64 + (a * d + b * c) * 2^32 + b * d
    $a = uint64(gmp_shiftr($multiplier, 32));
	
    $b = uint32($multiplier);
    $c = uint64(gmp_shiftr($multiplicand,32));
    $d = uint32($multiplicand);

	$ad = __emulu($a, $d);
	$bd = __emulu($b, $d);

	$adbc = uint64(gmp_add($ad,__emulu($b , $c)));
    
    if(gmp_cmp($adbc, $ad) < 0) {
        $adbc_carry = 1;
    } else {
        $adbc_carry = 0;
    }

	// multiplier * multiplicand = product_hi * 2^64 + product_lo
	$product_lo = uint64(gmp_add($bd,uint64(gmp_shiftl($adbc, 32))));

    if(gmp_cmp($product_lo, $bd) < 0) {
        $product_lo_carry = 1;
    } else {
        $product_lo_carry = 0;
    }

	$product_hi = uint64(gmp_add(uint64(gmp_add(uint64(gmp_add(__emulu($a , $c),uint64(gmp_shiftr($adbc, 32)))), uint64(gmp_shiftl($adbc_carry, 32)))), $product_lo_carry));

	return $product_lo;
}

function ui128_quotient_mod($lo, $hi) {
	// hi:lo * ceil(2**170/MOD) >> (64 + 64 + 42)
	$prod1 = gmp_init(0);
	__umul128($lo, gmp_init("0x604fa6a1c6346a87"), $prod1);
	$part1hi = gmp_init(0);
	$part1lo = __umul128($lo, gmp_init("0x2d351c6d04f8b"), $part1hi);
	$part2hi = gmp_init(0);
	$part2lo = __umul128($hi, gmp_init("0x604fa6a1c6346a87"), $part2hi);
	$sum1 = uint64(gmp_add($part1lo, $part2lo));

    $sum1carry = 0;

    if(gmp_cmp($sum1, $part1lo) < 0) {
	    $sum1carry = 1;
    } 

	$sum1 = uint64(gmp_add($sum1,$prod1));

    if(gmp_cmp($sum1, $prod1) < 0) {
	    $sum1carry = $sum1carry + 1;
    } 

    $sum1carry = gmp_init($sum1carry);

	$prod2 = uint64(gmp_add($part1hi, uint64(gmp_add($part2hi, $sum1carry))));
	$prod3hi = 0;
	$prod3lo = __umul128($hi, gmp_init("0x2d351c6d04f8b"), $prod3hi);
	$prod3lo = uint64(gmp_add($prod3lo, $prod2));
    if(gmp_cmp($prod3lo, $prod2) < 0) {
        $prod3hi = uint64(gmp_add($prod3hi, 1));
    } 
	
	return uint64(gmp_or(uint64(gmp_shiftr($prod3lo, 42)),uint64(gmp_shiftl($prod3hi, 22))));
}

function residue_mul($x, $y) { // there is a bug
    global $MOD;
    // * ceil(2**170/MOD) = 0x2d351 c6d04f8b|604fa6a1 c6346a87 for (p-1)*(p-1) max
	$hi = gmp_init(0);
	$lo = __umul128($x, $y, $hi);
	$quotient = ui128_quotient_mod($lo,$hi);
	return uint64(gmp_sub($lo, uint64(gmp_mul($quotient, $MOD))));
}

function residue_pow($x, $y) {
    if(gmp_cmp($y,0)==0) {
        return gmp_init(1);
    }
    
	$cur = $x;

    while(! (gmp_cmp(gmp_and($y,1),0)>0)) { //while (!(y & 1)) {
        $cur = residue_mul($cur, $cur);
        $y = uint64(gmp_shiftr($y, 1));
    }
	$res = $cur;

    while (gmp_cmp($y,0)!=0) { // while ((y >>= 1) != 0) {
        $y = uint64(gmp_shiftr($y, 1));

        $cur = residue_mul($cur, $cur);

        if(gmp_cmp(gmp_and($y, 1),0)>0) {
            $res = residue_mul($res, $cur);
        }
    }

    $y = uint64(gmp_shiftr($y, 1)); 
	return $res;
}

function inverse($u, $v) {
	$tmp = 0; // i64
	$xu = 1; // i64
    $xv = 0; // i64
	$v0 = $v;
	while (gmp_cmp($u,1)>0) {
		//ui64 d = v / u; ui64 remainder = v % u;
        $d = uint64(gmp_div($v,$u));
        $remainder = uint64(gmp_mod($v, $u));
		$tmp = int64($u); 
        $u = $remainder; 
        $v = uint64($tmp);
		$tmp = $xu; 
        $xu = int64(gmp_sub($xv,gmp_mul($d, $xu))); 
        $xv = int64($tmp);
	}
	//xu += (xu < 0 ? v0 : 0);
    if(gmp_cmp($xu,0)<0) {
        $xu = int64(gmp_add($xu,$v0));
    }
	return uint64($xu);
}

function residue_inv($x) { 
    global $MOD;
    return inverse($x, $MOD); 
}

// #define BAD 0xFFFFFFFFFFFFFFFFull
// #define NON_RESIDUE 43
$BAD = gmp_init("0xFFFFFFFFFFFFFFFF");
$NON_RESIDUE = gmp_init("43");

function residue_sqrt($what) {
    global $BAD, $MOD, $NON_RESIDUE, $debugEN;
    // if(!what) - if zero
	if (gmp_cmp($what,0)==0)
		return gmp_init(0);

	$g = $NON_RESIDUE;
    $z = gmp_init(0);
    $y = $z;
    $r = $z;
    $x = $z;
    $b = $z;
    $t = $z;

	$e = $z;
    $q = uint64(gmp_sub($MOD,1));
	//while (!(q & 1))
	//	e++, q >>= 1;

    while (! (gmp_cmp(gmp_and($q, 1),0)>0)) {
        $e = uint64(gmp_add($e, 1));
        $q = uint64(gmp_shiftr($q, 1));
    }

	$z = residue_pow($g, $q);
	$y = $z;
	$r = $e;
	$x = residue_pow($what, uint64(gmp_div(uint64(gmp_sub($q,1)),2))); // (q - 1) / 2
	$b = residue_mul(residue_mul($what, $x), $x);
	$x = residue_mul($what, $x);

	while (gmp_cmp($b,1)!=0) {
		$m = gmp_init(0);
        $b2 = $b;
		do {
			$m = uint64(gmp_add($m,1));
			$b2 = residue_mul($b2, $b2);
		} while (gmp_cmp($b2,1)!=0);
		if (gmp_cmp($m, $r)==0)
			return $BAD;
		$t = residue_pow($y, uint64(gmp_shiftl(1, uint64(gmp_sub(uint64(gmp_sub($r, $m)), 1)))));
		$y = residue_mul($t, $t);
		$r = $m;
		$x = residue_mul($x, $t);
		$b = residue_mul($b, $y);
	}
	if (gmp_cmp(residue_mul($x, $x),$what)!=0) {
		//printf("internal error in sqrt\n");
		return $BAD;
	}
	return $x;
}

/*
    typedef struct {
        ui64 u[2];
        ui64 v[2];
    } TDivisor;
*/

class TDivisor {
    public $u;
    public $v;
}

$f = Array(gmp_init(0), gmp_init("0x21840136C85381"),gmp_init("0x44197B83892AD0"),gmp_init("0x1400606322B3B04"),gmp_init("0x1400606322B3B04"),gmp_init("1"));

function find_divisor_v(&$d) { // int find_divisor_v(TDivisor* d), u[] input v[] output
    global $f, $BAD;
	// u | v^2 - f
	// u = u0 + u1*x + x^2
	// f%u = f0 + f1*x
	$v1 = gmp_init(0);
	$f2 = Array(); //ui64 f2[6];
	//int i, j;
    $i = 0;
    $j = 0;
	for ($i = 0; $i < 6; $i++) {
        $f2[$i] = $f[$i];
    }
		
	$u0 = $d->u[0];
	$u1 = $d->u[1];

	for($j=3; $j>=0; $j--) {
		$f2[$j] = residue_sub($f2[$j], residue_mul($u0, $f2[$j + 2]));
		$f2[$j + 1] = residue_sub($f2[$j + 1], residue_mul($u1, $f2[$j + 2]));
		$f2[$j + 2] = 0;
	}
	// v = v0 + v1*x
	// u | (v0^2 - f0) + (2*v0*v1 - f1)*x + v1^2*x^2 = u0*v1^2 + u1*v1^2*x + v1^2*x^2
	// v0^2 - f0 = u0*v1^2
	// 2*v0*v1 - f1 = u1*v1^2
	// v0^2 = f0 + u0*v1^2 = (f1 + u1*v1^2)^2 / (2*v1)^2
	// (f1^2) + 2*(f1*u1-2*f0) * v1^2 + (u1^2-4*u0) * v1^4 = 0
	// v1^2 = ((2*f0-f1*u1) +- 2*sqrt(-f0*f1*u1 + f0^2 + f1^2*u0))) / (u1^2-4*u0)
	$f0 = $f2[0];
	$f1 = $f2[1];
	$u0double = residue_add($u0, $u0);
	$coeff2 = residue_sub(residue_mul($u1, $u1), residue_add($u0double, $u0double));
	$coeff1 = residue_sub(residue_add($f0, $f0), residue_mul($f1, $u1));
	if (gmp_cmp($coeff2, 0) == 0) {
		if (gmp_cmp($coeff1, 0) == 0) {
			if (gmp_cmp($f1, 0) == 0) {
				// impossible
				//printf("bad f(), double root detected\n");
			}
			return 0;
		}
		$sqr = residue_mul(residue_mul($f1, $f1), residue_inv(residue_add($coeff1, $coeff1)));
		$v1 = residue_sqrt($sqr);
		if (gmp_cmp($v1, $BAD) == 0)
			return 0;
	} else {
		$d = residue_add(residue_mul($f0, $f0), residue_mul($f1, residue_sub(residue_mul($f1, $u0), residue_mul($f0, $u1))));
		$d = residue_sqrt($d);
		if (gmp_cmp($d, $BAD) == 0)
			return 0;
		$d = residue_add($d, $d);
		$inv = residue_inv($coeff2);
		$root = residue_mul(residue_add($coeff1, $d), $inv);
		$v1 = residue_sqrt($root);
		if (gmp_cmp($v1, $BAD) == 0) {
			$root = residue_mul(residue_sub($coeff1, $d), $inv);
			$v1 = residue_sqrt($root);
			if (gmp_cmp($v1, $BAD) == 0)
				return 0;
		}
	}
	$v0 = residue_mul(residue_add($f1, residue_mul($u1, residue_mul($v1, $v1))), residue_inv(residue_add($v1, $v1)));
	$d->u[0] = $u0;
	$d->u[1] = $u1;
	$d->v[0] = $v0;
	$d->v[1] = $v1;
	return 1;
}

// UNTESTED 
function polynomial_mul($adeg, $a, $bdeg, $b, $resultprevdeg, &$result) {  //static int polynomial_mul(int adeg, const ui64 a[], int bdeg, const ui64 b[], int resultprevdeg, ui64 result[])
    // $a[], $b[], $result[] is uint64
	if ($adeg < 0 || $bdeg < 0)
		return $resultprevdeg;
    //	int i, j;
    $i = 0;
    $j = 0;

	for ($i = $resultprevdeg + 1; $i <= $adeg + $bdeg; $i++)
		$result[$i] = gmp_init(0);
	$resultprevdeg = $i - 1;
	for ($i = 0; $i <= $adeg; $i++)
		for ($j = 0; $j <= $bdeg; $j++)
			$result[$i + $j] = residue_add($result[$i + $j], residue_mul($a[$i], $b[$j]));
	while ($resultprevdeg >= 0 && gmp_cmp($result[$resultprevdeg],0) == 0)
		--$resultprevdeg;
	return $resultprevdeg;
}

// UNTESTED
function polynomial_div_monic($adeg, &$a, $bdeg, $b, &$quotient) { //static int polynomial_div_monic(int adeg, ui64 a[], int bdeg, const ui64 b[], ui64* quotient)
    // $a[], $b[], $quotient are uint64
	//assert(bdeg >= 0);
	//assert(b[bdeg] == 1);
	//int i, j;
	//echo("bdeg: $bdeg\n");
	//debug_print_backtrace();
    $i = 0;
    $j = 0;
	for ($i = ($adeg - $bdeg); $i >= 0; $i--) {
		$q = $a[$i + $bdeg]; // uint64
		//if (gmp_cmp($quotient[0],0)>0) 
		if(gettype($quotient)=="array")
			$quotient[$i] = $q;
		for ($j = 0; $j < $bdeg; $j++)
			$a[$i + $j] = residue_sub($a[$i + $j], residue_mul($q, $b[$j]));
		$a[$i + $j] = gmp_init(0);
	}
	$i += $bdeg;
	while ($i >= 0 && gmp_cmp($a[$i], 0) == 0)
		$i--;
	return $i;
}

// UNTESTED
function polynomial_xgcd($adeg, $a, $bdeg, $b, &$pgcddeg, &$gcd, &$pmult1deg, &$mult1, &$pmult2deg, &$mult2) { //static void polynomial_xgcd(int adeg, const ui64 a[3], int bdeg, const ui64 b[3], int* pgcddeg, ui64 gcd[3], int* pmult1deg, ui64 mult1[3], int* pmult2deg, ui64 mult2[3])
    // $a[], $b[], $gcd[], $mult1[], $mult2[] are uint64
	$sdeg = -1;
	//ui64 s[3] = {0, 0, 0};
    $s = Array(gmp_init(0),gmp_init(0),gmp_init(0));
	$mult1deg = 0;
	$mult1[0] = gmp_init(1); 
    $mult1[1] = gmp_init(0);
    $mult1[2] = gmp_init(0);

	$tdeg = 0;
	//ui64 t[3] = {1, 0, 0};
    $t = Array(gmp_init(1),gmp_init(0),gmp_init(0));
	$mult2deg = -1;
	$mult2[0] = gmp_init(0);
    $mult2[1] = gmp_init(0); 
    $mult2[2] = gmp_init(0);
	$rdeg = $bdeg;
	//ui64 r[3] = {b[0], b[1], b[2]};
    $r = Array($b[0],$b[1],$b[2]);
	$gcddeg = $adeg;
	$gcd[0] = $a[0];
    $gcd[1] = $a[1];
    $gcd[2] = $a[2];
	// s*u1 + t*u2 = r
	// mult1*u1 + mult2*u2 = gcd
	while ($rdeg >= 0) {
		if ($rdeg > $gcddeg) {
			//unsigned tmp;
            $tmp = 0;
			//int tmpi;
            $tmpi = 0;

			$tmp = $rdeg; $rdeg = $gcddeg; $gcddeg = $tmp;
			$tmpi = $sdeg; $sdeg = $mult1deg; $mult1deg = $tmpi;
			$tmpi = $tdeg; $tdeg = $mult2deg; $mult2deg = $tmpi;
			//ui64 tmp2;
            $tmp2 = gmp_init(0);
			$tmp2 = $r[0]; $r[0] = $gcd[0]; $gcd[0] = $tmp2;
			$tmp2 = $r[1]; $r[1] = $gcd[1]; $gcd[1] = $tmp2;
			$tmp2 = $r[2]; $r[2] = $gcd[2]; $gcd[2] = $tmp2;
			$tmp2 = $s[0]; $s[0] = $mult1[0]; $mult1[0] = $tmp2;
			$tmp2 = $s[1]; $s[1] = $mult1[1]; $mult1[1] = $tmp2;
			$tmp2 = $s[2]; $s[2] = $mult1[2]; $mult1[2] = $tmp2;
			$tmp2 = $t[0]; $t[0] = $mult2[0]; $mult2[0] = $tmp2;
			$tmp2 = $t[1]; $t[1] = $mult2[1]; $mult2[1] = $tmp2;
			$tmp2 = $t[2]; $t[2] = $mult2[2]; $mult2[2] = $tmp2;
			continue;
		}
		//int delta = gcddeg - rdeg;
        $delta = $gcddeg - $rdeg;
		$mult = residue_mul($gcd[$gcddeg], residue_inv($r[$rdeg])); // ui64
		// quotient = mult * x**delta
		//assert(rdeg + delta < 3);
		for ($i = 0; $i <= $rdeg; $i++)
			$gcd[$i + $delta] = residue_sub($gcd[$i + $delta], residue_mul($mult, $r[$i]));
		while ($gcddeg >= 0 && gmp_cmp($gcd[$gcddeg], 0) == 0)
			$gcddeg--;
		//assert(sdeg + delta < 3);
		for ($i = 0; $i <= $sdeg; $i++)
			$mult1[$i + $delta] = residue_sub($mult1[$i + $delta], residue_mul($mult, $s[$i]));
		if ($mult1deg < $sdeg + $delta)
			$mult1deg = $sdeg + $delta;
		while ($mult1deg >= 0 && gmp_cmp($mult1[$mult1deg],0) == 0)
			$mult1deg--;
		//assert(tdeg + delta < 3);
		for ($i = 0; $i <= $tdeg; $i++)
			$mult2[$i + $delta] = residue_sub($mult2[$i + $delta], residue_mul($mult, $t[$i]));
		if ($mult2deg < $tdeg + $delta)
			$mult2deg = $tdeg + $delta;
		while ($mult2deg >= 0 && gmp_cmp($mult2[$mult2deg], 0) == 0)
			$mult2deg--;
	}
	// d1 = gcd, e1 = mult1, e2 = mult2
	$pgcddeg = $gcddeg;
	$pmult1deg = $mult1deg;
	$pmult2deg = $mult2deg;
}

// UNTESTED
function u2poly(&$src, &$polyu, &$polyv) { //static int u2poly(const TDivisor* src, ui64 polyu[3], ui64 polyv[2])
    // $polyu[] and $polyv[] are uint64
    global $BAD;

	if (gmp_cmp($src->u[1], $BAD) != 0) {
		$polyu[0] = $src->u[0];
		$polyu[1] = $src->u[1];
		$polyu[2] = 1;
		$polyv[0] = $src->v[0];
		$polyv[1] = $src->v[1];
		return 2;
	}
	if (gmp_cmp($src->u[0], $BAD) != 0) {
		$polyu[0] = $src->u[0];
		$polyu[1] = 1;
		$polyv[0] = $src->v[0];
		$polyv[1] = 0;
		return 1;
	}
	$polyu[0] = gmp_init(1);
	$polyv[0] = gmp_init(0);
	$polyv[1] = gmp_init(0);
	return 0;
}

// UNTESTED
function divisor_add($src1, $src2, &$dst) {// static void divisor_add(const TDivisor* src1, const TDivisor* src2, TDivisor* dst)
	//ui64 u1[3], u2[3], v1[2], v2[2];
    global $f, $BAD;
    $gmp0 = gmp_init(0);

    $u1 = Array($gmp0,$gmp0,$gmp0);
    $u2 = Array($gmp0,$gmp0,$gmp0);
    $v1 = Array($gmp0,$gmp0);
    $v2 = Array($gmp0,$gmp0);

	$u1deg = u2poly($src1, $u1, $v1); // int
	$u2deg = u2poly($src2, $u2, $v2); // int

	//debug_print_backtrace();
	

	// extended gcd: d1 = gcd(u1, u2) = e1*u1 + e2*u2
	//int d1deg, e1deg, e2deg;
    $d1deg = 0;
    $e1deg = 0;
    $e2deg = 0;
	//ui64 d1[3], e1[3], e2[3];
    $d1 = Array($gmp0,$gmp0,$gmp0);
    $e1 = Array($gmp0,$gmp0,$gmp0);
    $e2 = Array($gmp0,$gmp0,$gmp0);
	polynomial_xgcd($u1deg, $u1, $u2deg, $u2, $d1deg, $d1, $e1deg, $e1, $e2deg, $e2);
	
	//assert(e1deg <= 1);
	//assert(e2deg <= 1);
	// extended gcd again: d = gcd(d1, v1+v2) = c1*d1 + c2*(v1+v2)
	//ui64 b[3] = {residue_add(v1[0], v2[0]), residue_add(v1[1], v2[1]), 0};
    $b = Array(residue_add($v1[0], $v2[0]), residue_add($v1[1], $v2[1]), gmp_init(0));
	//int bdeg = (b[1] == 0 ? (b[0] == 0 ? -1 : 0) : 1); 
    $bdeg = 1;
    if(gmp_cmp($b[1], 0) == 0) {
        if(gmp_cmp($b[0], 0) == 0) {
            $bdeg = -1;
        } else {
            $bdeg = 0;
        }
    }
	//for($i=0;$i<3;$i++) {
	//	echo("b[$i]: ".strval($b[$i])."\n");
	//}

	//int ddeg, c1deg, c2deg;
    $ddeg = 0;
    $c1deg = 0;
    $c2deg = 0;
	//ui64 d[3], c1[3], c2[3];
    $d = Array($gmp0,$gmp0,$gmp0);
    $c1 = Array($gmp0,$gmp0,$gmp0);
    $c2 = Array($gmp0,$gmp0,$gmp0);
	polynomial_xgcd($d1deg, $d1, $bdeg, $b, $ddeg, $d, $c1deg, $c1, $c2deg, $c2);

	//echo("ddeg: $u1deg\n");
	//echo("c1deg: $c1deg\n");
	//echo("c2deg: $c2deg\n");

	/*
	for($i=0;$i<3;$i++) {
		echo("d[$i]: ".strval($d[$i])."\n");
	}
	for($i=0;$i<3;$i++) {
		echo("c1[$i]: ".strval($c1[$i])."\n");
	}
	for($i=0;$i<3;$i++) {
		echo("c2[$i]: ".strval($c2[$i])."\n");
	}*/

	//assert(c1deg <= 0);
	//assert(c2deg <= 1);
	//assert(ddeg >= 0);
	$dmult = residue_inv($d[$ddeg]); // ui64
	
	for ($i = 0; $i < $ddeg; $i++)
		$d[$i] = residue_mul($d[$i], $dmult);
	$d[$i] = 1;
	for ($i = 0; $i <= $c1deg; $i++)
		$c1[$i] = residue_mul($c1[$i], $dmult);
	for ($i = 0; $i <= $c2deg; $i++)
		$c2[$i] = residue_mul($c2[$i], $dmult);
	//ui64 u[5];
    $u = Array($gmp0,$gmp0,$gmp0,$gmp0,$gmp0);
	// echo("u1deg, u2deg: $u1deg $u2deg\n");
	// for($i=0;$i<3;$i++) {
	// 	echo("u1[$i]: ".strval($u1[$i])."\n");
	// }
	// for($i=0;$i<3;$i++) {
	// 	echo("u2[$i]: ".strval($u2[$i])."\n");
	// }
	
	$udeg = polynomial_mul($u1deg, $u1, $u2deg, $u2, -1, $u); // int
	// for($i=0;$i<5;$i++) {
	// 	echo("u[$i]: ".strval($u[$i])."\n");
	// }
	// echo("udeg: $udeg\n");
	// u is monic
	//ui64 v[7], tmp[7];
    $v = Array($gmp0,$gmp0,$gmp0,$gmp0,$gmp0,$gmp0,$gmp0);
    $tmp = Array($gmp0,$gmp0,$gmp0,$gmp0,$gmp0,$gmp0,$gmp0);
	//int vdeg, tmpdeg;
    $vdeg = 0;
    $tmpdeg = 0;
	// c1*(e1*u1*v2 + e2*u2*v1) + c2*(v1*v2 + f)
	// c1*(e1*u1*(v2-v1) + d1*v1) + c2*(v1*v2 + f)
	$v[0] = residue_sub($v2[0], $v1[0]);
	$v[1] = residue_sub($v2[1], $v1[1]);
	$tmpdeg = polynomial_mul($e1deg, $e1, 1, $v, -1, $tmp);
	$vdeg = polynomial_mul($u1deg, $u1, $tmpdeg, $tmp, -1, $v);
	$vdeg = polynomial_mul($d1deg, $d1, 1, $v1, $vdeg, $v);
	for ($i = 0; $i <= $vdeg; $i++)
		$v[$i] = residue_mul($v[$i], $c1[0]);
	//memcpy(tmp, f, 6 * sizeof(f[0]));
    for($iii = 0; $iii < 6; $iii++) {
        $tmp[$iii] = $f[$iii];
    }

	$tmpdeg = 5;
	$tmpdeg = polynomial_mul(1, $v1, 1, $v2, $tmpdeg, $tmp);
	$vdeg = polynomial_mul($c2deg, $c2, $tmpdeg, $tmp, $vdeg, $v);
	//echo("udeg, ddeg: $udeg $ddeg\n");
	if ($ddeg > 0) {
		//assert(udeg >= 2*ddeg);
		//ui64 udiv[5];
        $udiv = Array($gmp0,$gmp0,$gmp0,$gmp0,$gmp0);
		polynomial_div_monic($udeg, $u, $ddeg, $d, $udiv); $udeg -= $ddeg;
		polynomial_div_monic($udeg, $udiv, $ddeg, $d, $u); $udeg -= $ddeg;
		if ($vdeg >= 0) {
			//assert(vdeg >= ddeg);
			polynomial_div_monic($vdeg, $v, $ddeg, $d, $udiv); $vdeg -= $ddeg;
			//memcpy(v, udiv, (vdeg + 1) * sizeof(v[0])); what the?
            for($iii = 0; $iii < ($vdeg + 1); $iii++) {
                $v[$iii] = $udiv[$iii];
            }
		}
	}
    $empty = gmp_init(0); // NULL placeholder
	$vdeg = polynomial_div_monic($vdeg, $v, $udeg, $u, $empty);
	while ($udeg > 2) {
		//assert(udeg <= 4);
		//assert(vdeg <= 3);
		// u' = monic((f-v^2)/u), v'=-v mod u'
		$tmpdeg = polynomial_mul($vdeg, $v, $vdeg, $v, -1, $tmp);
		for ($i = 0; $i <= $tmpdeg && $i <= 5; $i++)
			$tmp[$i] = residue_sub($f[$i], $tmp[$i]);
		for (; $i <= $tmpdeg; $i++)
			$tmp[$i] = residue_sub(0, $tmp[$i]);
		for (; $i <= 5; $i++)
			$tmp[$i] = $f[$i];
		$tmpdeg = $i - 1;
		$udiv = Array($gmp0,$gmp0,$gmp0,$gmp0,$gmp0);
		polynomial_div_monic($tmpdeg, $tmp, $udeg, $u, $udiv);
		$udeg = $tmpdeg - $udeg;
		$mult = residue_inv($udiv[$udeg]); // uint64
		for ($i = 0; $i < $udeg; $i++)
			$u[$i] = residue_mul($udiv[$i], $mult);
		$u[$i] = 1;
		for ($i = 0; $i <= $vdeg; $i++)
			$v[$i] = residue_sub(0, $v[$i]);
        $empty = 0;
		$vdeg = polynomial_div_monic($vdeg, $v, $udeg, $u, $empty);
	}
	if ($udeg == 2) {
		$dst->u[0] = $u[0];
		$dst->u[1] = $u[1];
		//$dst->v[0] = (vdeg >= 0 ? v[0] : 0);
        if($vdeg>=0) {
            $dst->v[0] = $v[0];
        } else {
            $dst->v[0] = 0;
        }
		//$dst->v[1] = (vdeg >= 1 ? v[1] : 0);
        if($vdeg>=1) {
            $dst->v[1] = $v[1];
        } else {
            $dst->v[1] = 0;
        }
	} else if ($udeg == 1) {
		$dst->u[0] = $u[0];
		$dst->u[1] = $BAD;
		//$dst->v[0] = (vdeg >= 0 ? v[0] : 0);
        if($vdeg>=0) {
            $dst->v[0] = $v[0];
        } else {
            $dst->v[0] = 0;
        }
		$dst->v[1] = $BAD;
	} else {
		//assert(udeg == 0);
		$dst->u[0] = $BAD;
		$dst->u[1] = $BAD;
		$dst->v[0] = $BAD;
		$dst->v[1] = $BAD;
	}
}

// UNTESTED
// #define divisor_double(src, dst) divisor_add(src, src, dst)
function divisor_double($src, &$dst) {
    return divisor_add($src, $src, $dst);
}

// UNTESTED
function divisor_mul(&$src, $mult, &$dst) { //static void divisor_mul(const TDivisor* src, ui64 mult, TDivisor* dst)
    // mult is uint64
    global $BAD;

	if (gmp_cmp($mult, 0) == 0) {
		$dst->u[0] = $BAD;
		$dst->u[1] = $BAD;
		$dst->v[0] = $BAD;
		$dst->v[1] = $BAD;
		return;
	}
	//TDivisor cur = *src; this defines cur and copies src to cur
    $cur = new TDivisor();
    $cur->u = Array($src->u[0],$src->u[1]);
    $cur->v = Array($src->v[0],$src->v[1]);

	while (!(gmp_cmp(gmp_and($mult, 1),0)>0)) {
		divisor_double($cur, $cur);
        $mult = uint64(gmp_shiftr($mult,1));
	}
	// *dst = cur; copy cur to dst
    $dst->u = Array($cur->u[0],$cur->u[1]);
    $dst->v = Array($cur->v[0],$cur->v[1]);

    $mult = uint64(gmp_shiftr($mult,1));
	while(gmp_cmp($mult, 0) != 0) {// while ((mult >>= 1) != 0) {
		divisor_double($cur, $cur);
		if (gmp_cmp(gmp_and($mult, 1), 0) > 0)
			divisor_add($dst, $cur, $dst);
        $mult = uint64(gmp_shiftr($mult,1));
	}
}

// UNTESTED
function divisor_mul128($src, $mult_lo, $mult_hi, &$dst) { // static void divisor_mul128(const TDivisor* src, ui64 mult_lo, ui64 mult_hi, TDivisor* dst)
    // mult_lo and mult_hi are uint64
    global $BAD;

	if (gmp_cmp($mult_lo,0) == 0 && gmp_cmp($mult_hi, 0) == 0) {
		$dst->u[0] = $BAD;
		$dst->u[1] = $BAD;
		$dst->v[0] = $BAD;
		$dst->v[1] = $BAD;
		return;
	}
	//TDivisor cur = *src;
    $cur = new TDivisor();
	
    $cur->u = Array($src->u[0],$src->u[1]);
    $cur->v = Array($src->v[0],$src->v[1]);


	while (!(gmp_cmp(gmp_and($mult_lo, 1),0)>0)) {
		divisor_double($cur, $cur);
		
		$mult_lo = uint64(gmp_shiftr($mult_lo, 1));
		if (gmp_cmp(gmp_and($mult_hi, 1), 0) > 0)
			$mult_lo = uint64(gmp_or($mult_lo,uint64(gmp_shiftl(1,63))));
		$mult_hi = uint64(gmp_shiftr($mult_hi, 1));
	}
    //	*dst = cur;
    $dst->u = Array($cur->u[0],$cur->u[1]);
    $dst->v = Array($cur->v[0],$cur->v[1]);

	

	for (;;) {
		$mult_lo = uint64(gmp_shiftr($mult_lo, 1));
		if (gmp_cmp(gmp_and($mult_hi, 1), 0) > 0)
            $mult_lo = uint64(gmp_or($mult_lo,uint64(gmp_shiftl(1,63))));
            $mult_hi = uint64(gmp_shiftr($mult_hi, 1));
		if (gmp_cmp($mult_lo, 0) == 0 && gmp_cmp($mult_hi, 0) == 0)
			break;
		divisor_double($cur, $cur);
		//echo("\nu0: " . gmp_strval($cur->u[0]));
		//echo("\nu1: " . gmp_strval($cur->u[1]));
		//echo("\nv0: " . gmp_strval($cur->v[0]));
		//echo("\nv1: " . gmp_strval($cur->v[1]) . "\n");
		
		if (gmp_cmp(gmp_and($mult_lo, 1), 0) > 0) {
			divisor_add($dst, $cur, $dst);
		}
	}
}

function rol($x, $shift) {
	//assert(shift > 0 && shift < 32);
    if(gettype($shift)=="object") { 
        $shift = gmp_intval($shift);
    }
	return uint32(gmp_or(gmp_shiftl($x, $shift), gmp_shiftr($x, (32 - $shift))));
}

// array fun yay D: !!
// UNTESTED
function sha1_single_block($input, &$output) { //static void sha1_single_block(unsigned char input[64], unsigned char output[20])
	//unsigned a, b, c, d, e;
	$a = gmp_init("0x67452301");
	$b = gmp_init("0xEFCDAB89");
	$c = gmp_init("0x98BADCFE");
	$d = gmp_init("0x10325476");
	$e = gmp_init("0xC3D2E1F0");
	//unsigned w[80];
    $w = Array();
	
	for ($i = 0; $i < 16; $i++) { // ok
		$w[$i] = uint32(gmp_or(gmp_or(gmp_or(gmp_shiftl($input[4*$i], 24), gmp_shiftl($input[4*$i+1], 16)), gmp_shiftl($input[4*$i+2], 8)), $input[4*$i+3]));	
	}
	for ($i = 16; $i < 80; $i++) { // ok
		$w[$i] = rol(gmp_xor(gmp_xor(gmp_xor($w[$i - 3], $w[$i - 8]), $w[$i - 14]), $w[$i - 16]), 1);
		//echo($w[$i] . " ");
	}
	
	for ($i = 0; $i < 20; $i++) { // ?
        // for ~b
		//echo("A STEP1: $a\n");
		// bitwse not on b
        $b_inv = gmp_intval($b);
		$b_inv = ~$b_inv;
		$b_inv = $b_inv & 0xFFFFFFFF;
		$b_inv = gmp_init($b_inv);
		//$tmp = rol(a, tmp) + ((b & c) | (~b & d)) + e + w[i] + 0x5A827999;
		$tmp = uint32(gmp_or(gmp_and($b,$c),gmp_and($b_inv,$d)));
		//echo("A STEP2: $tmp\n");
		$tmp = uint32(gmp_add($tmp, rol($a, 5)));
		//echo("A STEP3: $tmp\n");
		$tmp = uint32(gmp_add($tmp,$e));
		//echo("A STEP4: $tmp\n");
		$tmp = uint32(gmp_add($tmp,$w[$i]));
		//echo("A STEP5: $tmp\n");
		$tmp = uint32(gmp_add($tmp,gmp_init("0x5A827999")));
		//$tmp = uint32(gmp_add(gmp_add(gmp_add(rol($a, 5),gmp_or(gmp_and($b,$c),gmp_and($b_inv, $d))),$e),$w[$i]),gmp_init("0x5A827999"));
        //$tmp = uint32(gmp_add(uint32(gmp_add(uint32(gmp_add(gmp_add(rol($a, 5),gmp_or(gmp_and($b, $c), gmp_and($b_inv, $d))),$e)),$w[$i])),gmp_init("0x5A827999")) );
		$e = $d;
		$d = $c;
		$c = uint32(rol($b, 30));
		$b = $a;
		$a = $tmp;
		
	}
	// echo("\n");
	// echo("A: " . gmp_strval($a) . "\n");
	// echo("B: " . gmp_strval($b) . "\n");
	// echo("C: " . gmp_strval($c) . "\n");
	// echo("D: " . gmp_strval($d) . "\n");
	// echo("E: " . gmp_strval($e) . "\n");
	
	for ($i = 20; $i < 40; $i++) {
		//unsigned tmp = rol(a, 5) + (b ^ c ^ d) + e + w[i] + 0x6ED9EBA1;
        $tmp = uint32(gmp_add(gmp_add(gmp_add(gmp_add(rol($a, 5),gmp_xor(gmp_xor($b,$c),$d)),$e),$w[$i]),gmp_init("0x6ED9EBA1")));
		$e = $d;
		$d = $c;
		$c = rol($b, 30);
		$b = $a;
		$a = $tmp;
	}
	for ($i = 40; $i < 60; $i++) {
		//unsigned tmp = rol(a, 5) + ((b & c) | (b & d) | (c & d)) + e + w[i] + 0x8F1BBCDC;
        $tmp = uint32(gmp_add(gmp_add(gmp_add(gmp_add(rol($a,5),gmp_or(gmp_or(gmp_and($b, $c), gmp_and($b, $d)), gmp_and($c, $d))),$e),$w[$i]),gmp_init("0x8F1BBCDC")));
		$e = $d;
		$d = $c;
		$c = rol($b, 30);
		$b = $a;
		$a = $tmp;
	}
	for ($i = 60; $i < 80; $i++) {
		//unsigned tmp = rol(a, 5) + (b ^ c ^ d) + e + w[i] + 0xCA62C1D6;
        $tmp = uint32(gmp_add(gmp_add(gmp_add(gmp_add(rol($a,5), gmp_xor(gmp_xor($b, $c),$d)), $e), $w[$i]), gmp_init("0xCA62C1D6")));
		$e = $d;
		$d = $c;
		$c = rol($b, 30);
		$b = $a;
		$a = $tmp;
	}
	$a = uint32(gmp_add($a, gmp_init("0x67452301")));
	$b = uint32(gmp_add($b, gmp_init("0xEFCDAB89")));
	$c = uint32(gmp_add($c, gmp_init("0x98BADCFE")));
	$d = uint32(gmp_add($d, gmp_init("0x10325476")));
	$e = uint32(gmp_add($e, gmp_init("0xC3D2E1F0")));
	
	$output[0] = uint8(gmp_shiftr($a, 24)); $output[1] = uint8(gmp_shiftr($a, 16)); $output[2] = uint8(gmp_shiftr($a, 8)); $output[3] = uint8($a);
	$output[4] = uint8(gmp_shiftr($b, 24)); $output[5] = uint8(gmp_shiftr($b, 16)); $output[6] = uint8(gmp_shiftr($b, 8)); $output[7] = uint8($b);
	$output[8] = uint8(gmp_shiftr($c, 24)); $output[9] = uint8(gmp_shiftr($c, 16)); $output[10] = uint8(gmp_shiftr($c, 8)); $output[11] = uint8($c);
	$output[12] = uint8(gmp_shiftr($d, 24)); $output[13] = uint8(gmp_shiftr($d, 16)); $output[14] = uint8(gmp_shiftr($d, 8)); $output[15] = uint8($d);
	$output[16] = uint8(gmp_shiftr($e, 24)); $output[17] = uint8(gmp_shiftr($e, 16)); $output[18] = uint8(gmp_shiftr($e, 8)); $output[19] = uint8($e);
}

// UNTESTED
function Mix(&$buffer, $bufSize, &$key, $keySize) { //static void Mix(unsigned char* buffer, size_t bufSize, const unsigned char* key, size_t keySize)
	//unsigned char sha1_input[64];
	$sha1_input = Array();
	//unsigned char sha1_result[20];
	$sha1_result = Array();
	$half = floor($bufSize / 2);
	//assert(half <= sizeof(sha1_result) && half + keySize <= sizeof(sha1_input) - 9);
	$external_counter = 0;
	for ($external_counter = 0; $external_counter < 4; $external_counter++) {
		//memset(sha1_input, 0, sizeof(sha1_input));
		for($iii=0;$iii<64;$iii++) {
			$sha1_input[$iii] = gmp_init(0);;
		}
		//memcpy(sha1_input, buffer + half, half);
		for($iii=0;$iii<$half;$iii++) {
			$sha1_input[$iii] = $buffer[$iii+$half];
		}
		//memcpy(sha1_input + half, key, keySize);
		for($iii=0;$iii<$keySize;$iii++) {
			$sha1_input[$iii + $half] = $key[$iii];
		}
		//sha1_input[half + keySize] = 0x80;
		$sha1_input[$half + $keySize] = gmp_init(0x80);
		//sha1_input[sizeof(sha1_input) - 1] = (half + keySize) * 8;
		$sha1_input[64 - 1] = uint8(gmp_mul(gmp_init(($half + $keySize)), 8));
		$sha1_input[64 - 2] = uint8(gmp_div(gmp_mul(gmp_init($half + $keySize), 8), 0x100));


		sha1_single_block($sha1_input, $sha1_result);
		
		for ($i = $half & ~3; $i < $half; $i++)
			$sha1_result[$i] = $sha1_result[$i + 4 - ($half & 3)];
		for ($i = 0; $i < $half; $i++) {
			$tmp = $buffer[$i + $half];
			$buffer[$i + $half] = gmp_xor($buffer[$i], $sha1_result[$i]);
			$buffer[$i] = $tmp;
		}
	}
}

// UNTESTED
function Unmix(&$buffer, $bufSize, &$key, $keySize) { // static void Unmix(unsigned char* buffer, size_t bufSize, const unsigned char* key, size_t keySize)
	//unsigned char sha1_input[64];
	$sha1_input = Array();
	//unsigned char sha1_result[20];
	$sha1_result = Array();
	$half = floor($bufSize / 2);

	//assert(half <= sizeof(sha1_result) && half + keySize <= sizeof(sha1_input) - 9);
	
	for ($external_counter = 0; $external_counter < 4; $external_counter++) {
		//memset(sha1_input, 0, sizeof(sha1_input));
		for($iii=0;$iii<64;$iii++) {
			$sha1_input[$iii] = gmp_init(0);
		}
		//memcpy(sha1_input, buffer, half);
		for($iii=0;$iii<$half;$iii++) {
			$sha1_input[$iii] = $buffer[$iii];
		}
		//memcpy(sha1_input + half, key, keySize);
		for($iii=0;$iii<$keySize;$iii++) {
			$sha1_input[$iii+$half] = $key[$iii];
		}
		$sha1_input[$half + $keySize] = gmp_init(0x80);
		$sha1_input[64 - 1] = uint8(gmp_mul(gmp_init(($half + $keySize)), 8));
		$sha1_input[64 - 2] = uint8(gmp_div(gmp_mul(gmp_init($half + $keySize), 8), 0x100));


		sha1_single_block($sha1_input, $sha1_result);


		for ($i = $half & ~3; $i < $half; $i++) {
			$sha1_result[$i] = $sha1_result[$i + 4 - ($half & 3)];
		}
		for ($i = 0; $i < $half; $i++) {
			$tmp = $buffer[$i];
			$buffer[$i] = gmp_xor($buffer[$i + $half], $sha1_result[$i]);
			$buffer[$i + $half] = $tmp;
		}
	}
}

$ERR_TOO_SHORT = 1;
$ERR_TOO_LARGE = 2;
$ERR_INVALID_CHARACTER = 3;
$ERR_INVALID_CHECK_DIGIT = 4;
$ERR_UNKNOWN_VERSION = 5;
$ERR_UNLUCKY = 6;

// UNTESTED
function generate($installation_id_str, &$confirmation_id) { // static int generate(const CHARTYPE* installation_id_str, CHARTYPE confirmation_id[49])
	global $ERR_TOO_SHORT, $ERR_TOO_LARGE, $ERR_INVALID_CHARACTER, $ERR_INVALID_CHECK_DIGIT, $ERR_UNKNOWN_VERSION, $ERR_UNLUCKY, $NON_RESIDUE, $MOD, $BAD;
	$installation_id = Array(); //unsigned char installation_id[19]; // 10**45 < 256**19
	$installation_id_len = 0; //size_t installation_id_len = 0;
	$p = $installation_id_str;
	$count = 0;
	$totalCount = 0;
	$check = 0; // unsigned
	$plen = strlen($p);
	for($pi=0;$pi<$plen;$pi++) { //for (; *p; p++) {
		$pp = $p[$pi];
		if ($pp == ' ' || $pp == '-')
			continue;
		/*
		int d = *p - '0';
		if (d < 0 || d > 9)
			return ERR_INVALID_CHARACTER;
		*/
		if(!is_numeric($pp)) {
			return $ERR_INVALID_CHARACTER;
		}
		$d = intval($pp);
		if ($count == 5 || ord($p[$pi+1]) == 0) {
			if (!$count)
				return ($totalCount == 45) ? $ERR_TOO_LARGE : $ERR_TOO_SHORT;
			if ($d != $check % 7)
				return ($count < 5) ? $ERR_TOO_SHORT : $ERR_INVALID_CHECK_DIGIT;
			$check = 0;
			$count = 0;
			continue;
		}
		$check += ($count % 2 ? $d * 2 : $d);
		$count++;
		$totalCount++;
		if ($totalCount > 45)
			return $ERR_TOO_LARGE;
		$carry = uint8(gmp_init($d)); // unsigned char
		for ($i = 0; $i < $installation_id_len; $i++) {
			$x = uint32(gmp_add(gmp_mul(gmp_init($installation_id[$i]), 10), $carry)); // unsigned
			$installation_id[$i] = gmp_intval(gmp_and($x, 0xFF));
			$carry = gmp_intval(uint8(gmp_shiftr($x, 8)));
		}
		if (gmp_cmp($carry, 0) > 0) {
			//assert(installation_id_len < sizeof(installation_id));
			$installation_id[$installation_id_len++] = gmp_intval($carry);
		}
	}
	if ($totalCount != 41 && $totalCount < 45)
		return $ERR_TOO_SHORT;
	for (; $installation_id_len < 19; $installation_id_len++)
		$installation_id[$installation_id_len] = 0;
	//static const unsigned char iid_key[4] = {0x6A, 0xC8, 0x5E, 0xD4};
	$iid_key = Array(gmp_init(0x6A), gmp_init(0xC8), gmp_init(0x5E), gmp_init(0xD4));
	// convert instalation id to gmp
	for($iii=0;$iii<19;$iii++) {
		$installation_id[$iii] = uint8(gmp_init($installation_id[$iii]));
	}

	Unmix($installation_id, $totalCount == 41 ? 17 : 19, $iid_key, 4);
	// convert it back to num (no need)
	//for($iii=0;$iii<19;$iii++) {
	//	$installation_id[$iii] = gmp_intval($installation_id[$iii]);
	//}


	if (gmp_cmp($installation_id[18], 0x10) >= 0)
		return $ERR_UNKNOWN_VERSION;

	// #pragma pack(push, 1)
	/*
	struct {
		ui64 HardwareID;
		ui64 ProductIDLow;
		unsigned char ProductIDHigh;
		unsigned short KeySHA1;
	} parsed;
	*/
	$parsed = new stdClass(); 
	$parsed->HardwareID = gmp_init(0); // uint64
	$parsed->ProductIDLow = gmp_init(0); // uint64
	$parsed->ProductIDHigh = gmp_init(0); // uint8
	$parsed->KeySHA1 = gmp_init(0); // uint16
	// #pragma pack(pop)
	//memcpy(&parsed, installation_id, sizeof(parsed));
	// splitted into multiple lines to not make very long lines
	$parsed->HardwareID = gmp_or(gmp_or(gmp_or($installation_id[0],gmp_shiftl($installation_id[1],8)),gmp_shiftl($installation_id[2],16)),gmp_shiftl($installation_id[3],24));
	$parsed->HardwareID = gmp_or($parsed->HardwareID,gmp_or(gmp_or(gmp_or(gmp_shiftl($installation_id[4],32),gmp_shiftl($installation_id[5],40)),gmp_shiftl($installation_id[6],48)),gmp_shiftl($installation_id[7],56)));

	$parsed->ProductIDLow = gmp_or(gmp_or(gmp_or($installation_id[8],gmp_shiftl($installation_id[9],8)),gmp_shiftl($installation_id[10],16)),gmp_shiftl($installation_id[11],24));
	$parsed->ProductIDLow = gmp_or($parsed->ProductIDLow,gmp_or(gmp_or(gmp_or(gmp_shiftl($installation_id[12],32),gmp_shiftl($installation_id[13],40)),gmp_shiftl($installation_id[14],48)),gmp_shiftl($installation_id[15],56)));

	$parsed->ProductIDHigh = $installation_id[16];

	$parsed->KeySHA1 = gmp_or($installation_id[17], gmp_shiftl($installation_id[18],8));
	// end of memcpy


	$productId1 = uint32(gmp_and($parsed->ProductIDLow, gmp_init((1 << 17) - 1))); // unsigned
	$productId2 = uint32(gmp_and(gmp_shiftr($parsed->ProductIDLow, 17), gmp_init((1 << 10) - 1))); // unsigned
	$productId3 = uint32(gmp_and(gmp_shiftr($parsed->ProductIDLow, 27), gmp_init((1 << 25) - 1))); // unsigned
	$version = uint32(gmp_and(gmp_shiftr($parsed->ProductIDLow, 52), 7)); // unsigned
	$productId4 = uint32(gmp_or(gmp_shiftr($parsed->ProductIDLow, 55), gmp_shiftl($parsed->ProductIDHigh, 9))); // unsigned
	if (gmp_cmp($version, gmp_init($totalCount == 41 ? 4 : 5))!=0)
		return $ERR_UNKNOWN_VERSION;
	//printf("Product ID: %05u-%03u-%07u-%05u\n", productId1, productId2, productId3, productId4);


	$keybuf = Array(); //unsigned char keybuf[16]; // unsigned char
	//memcpy(keybuf, &parsed.HardwareID, 8);
	// parsed->HardwareID is basically $installation_id[0 .. 7]
	for($iii=0;$iii<8;$iii++) {
		$keybuf[$iii] = $installation_id[$iii];
	}

	$productIdMixed = gmp_or(gmp_or(gmp_or(uint64(gmp_shiftl($productId1, 41)), uint64(gmp_shiftl($productId2, 58))), uint64(gmp_shiftl($productId3, 17))), uint64($productId4)); // uint64
	//memcpy(keybuf + 8, &productIdMixed, 8);
	// 0x1122334455667788 -> array[0] = 0x88, array[1] = 0x77 and on...
	
	for($iii=0;$iii<8;$iii++) {
		$keybuf[$iii+8] = uint8(gmp_and(gmp_shiftr($productIdMixed, $iii*8), gmp_init(0xFF)));
	}

	
	//TDivisor d;
	$d = new TDivisor();
    $d->u = Array(gmp_init(0),gmp_init(0));
    $d->v = Array(gmp_init(0),gmp_init(0));

	//$attempt - unsigned char
	for ($attempt = 0; $attempt <= 0x80; $attempt++) {
		//echo("attempt: $attempt\n");
		/*union {
			unsigned char buffer[14];
			struct {
				ui64 lo;
				ui64 hi;
			};
		} u;
		union makes char and the both ui64 be at the same memory
		u.buffer[0] = 0x64 makes u.lo = 0x64
		*/
		$u = new stdClass(); 
		$u->lo = gmp_init(0); 
		$u->hi = gmp_init(0);
		$u->buffer = Array();
		for($iii=0;$iii<14;$iii++) {
			$u->buffer[$iii] = gmp_init(0);
		}
		$u->buffer[7] = uint8(gmp_init($attempt));

		Mix($u->buffer, 14, $keybuf, 16);

		// copy buffer to lo and hi
		$u->lo = gmp_or(gmp_or(gmp_or($u->buffer[0],gmp_shiftl($u->buffer[1],8)),gmp_shiftl($u->buffer[2],16)),gmp_shiftl($u->buffer[3],24));
		$u->lo = gmp_or($u->lo,gmp_or(gmp_or(gmp_or(gmp_shiftl($u->buffer[4],32),gmp_shiftl($u->buffer[5],40)),gmp_shiftl($u->buffer[6],48)),gmp_shiftl($u->buffer[7],56)));

		$u->hi = gmp_or(gmp_or(gmp_or($u->buffer[8],gmp_shiftl($u->buffer[9],8)),gmp_shiftl($u->buffer[10],16)),gmp_shiftl($u->buffer[11],24));
		$u->hi = gmp_or($u->hi,gmp_or(gmp_shiftl($u->buffer[12],32),gmp_shiftl($u->buffer[13],40)));

		$x2 = ui128_quotient_mod($u->lo, $u->hi); //ui64
		$x1 = uint64(gmp_sub($u->lo, uint64(gmp_mul($x2, $MOD)))); // ui64
		$x2 = uint64(gmp_add($x2, 1));//x2++;
		$d->u[0] = residue_sub(residue_mul($x1, $x1), residue_mul($NON_RESIDUE, residue_mul($x2, $x2)));
		$d->u[1] = residue_add($x1, $x1);
		if (find_divisor_v($d))
			break;
	}
	if ($attempt > 0x80)
		return $ERR_UNLUCKY;

		
	divisor_mul128($d, gmp_init("0x04e21b9d10f127c1"), gmp_init("0x40da7c36d44c"), $d);


	/*union {
		struct {
			ui64 encoded_lo, encoded_hi;
		};
		struct {
			uint32_t encoded[4];
		};
	} e;
	*/
	$e = new stdClass(); 
	$e->encoded_lo = gmp_init(0);
	$e->encoded_hi = gmp_init(0);
	$e->encoded = Array(gmp_init(0),gmp_init(0),gmp_init(0),gmp_init(0));

	if ($d->u[0] == $BAD) {
		// we can not get the zero divisor, actually...
		$e->encoded_lo = __umul128(uint64(gmp_add($MOD, gmp_init(2))), $MOD, $e->encoded_hi);
	} else if ($d->u[1] == $BAD) {
		// O(1/MOD) chance
		//encoded = (unsigned __int128)(MOD + 1) * d.u[0] + MOD; // * MOD + d.u[0] is fine too
		$e->encoded_lo = __umul128(uint64(gmp_add($MOD, 1)), $d->u[0], $e->encoded_hi);
		$e->encoded_lo = uint64(gmp_add($e->encoded_lo, $MOD));
		if(gmp_cmp($e->encoded_lo, $MOD)<0) { // e.encoded_hi += (e.encoded_lo < MOD);
			$e->encoded_hi = uint64(gmp_add($e->encoded_hi, 1));
		}
		
	} else {
		//x1 = (d.u[1] % 2 ? d.u[1] + MOD : d.u[1]) / 2; // ui64
		if(gmp_cmp(gmp_mod($d->u[1],gmp_init(2)),0)>0) {
			$x1 = uint64(gmp_div(uint64(gmp_add($d->u[1], $MOD)),2));
		} else {
			$x1 = uint64(gmp_div($d->u[1],2));
		}
		$x2sqr = residue_sub(residue_mul($x1, $x1), $d->u[0]); // ui64
		$x2 = residue_sqrt($x2sqr); // ui64
		if (gmp_cmp($x2, $BAD) == 0) {
			$x2 = residue_sqrt(residue_mul($x2sqr, residue_inv($NON_RESIDUE)));
			//assert(x2 != BAD);
			$e->encoded_lo = __umul128(uint64(gmp_add($MOD, 1)), uint64(gmp_add($MOD, $x2)), $e->encoded_hi);
			$e->encoded_lo = uint64(gmp_add($e->encoded_lo, $x1));
			if(gmp_cmp($e->encoded_lo, $x1)<0) { // e.encoded_hi += (e.encoded_lo < x1);
				$e->encoded_hi = uint64(gmp_add($e->encoded_hi,1));
			}
			
		} else {
			// points (-x1+x2, v(-x1+x2)) and (-x1-x2, v(-x1-x2))
			// all ui64
			$x1a = residue_sub($x1, $x2);
			$y1 = residue_sub($d->v[0], residue_mul($d->v[1], $x1a));
			$x2a = residue_add($x1, $x2);
			$y2 = residue_sub($d->v[0], residue_mul($d->v[1], $x2a));
			if (gmp_cmp($x1a, $x2a)>0) {
				$tmp = $x1a; //ui64
				$x1a = $x2a;
				$x2a = $tmp;
			}
			//if ((y1 ^ y2) & 1) {
			if(gmp_cmp(gmp_and(gmp_xor($y1, $y2), 1), 0)>0) {
				$tmp = $x1a;
				$x1a = $x2a;
				$x2a = $tmp;
			}
			$e->encoded_lo = __umul128(uint64(gmp_add($MOD, 1)), $x1a, $e->encoded_hi);
			$e->encoded_lo = uint64(gmp_add($e->encoded_lo, $x2a));
			//e.encoded_hi += (e.encoded_lo < x2a);
			if(gmp_cmp($e->encoded_lo, $x2a) < 0) {
				$e->encoded_hi = uint64(gmp_add($e->encoded_hi, 1));
			}
		}
	}
	$decimal = Array();//unsigned char decimal[35];
	// convert e->encoeded_lo and e->encoded_hi to e->encoded
	$e->encoded[0] = gmp_and($e->encoded_lo,gmp_init("0xFFFFFFFF"));
	$e->encoded[1] = gmp_and(gmp_shiftr($e->encoded_lo,32),gmp_init("0xFFFFFFFF"));
	$e->encoded[2] = gmp_and($e->encoded_hi,gmp_init("0xFFFFFFFF"));
	$e->encoded[3] = gmp_and(gmp_shiftr($e->encoded_hi,32),gmp_init("0xFFFFFFFF"));

	for ($i = 0; $i < 35; $i++) {
		$c = uint32(gmp_mod($e->encoded[3], 10)); // unsigned
		$e->encoded[3] = uint32(gmp_div($e->encoded[3],10));
		//unsigned c2 = ((ui64)c << 32 | e.encoded[2]) % 10;
		$c2 = uint32(gmp_mod(gmp_or(uint64(gmp_shiftl($c,32)), $e->encoded[2]),10));
		//e.encoded[2] = ((ui64)c << 32 | e.encoded[2]) / 10;
		$e->encoded[2] = uint32(gmp_div(gmp_or(uint64(gmp_shiftl($c,32)), $e->encoded[2]),10));
		//unsigned c3 = ((ui64)c2 << 32 | e.encoded[1]) % 10;
		$c3 = uint32(gmp_mod(gmp_or(uint64(gmp_shiftl($c2,32)), $e->encoded[1]),10));
		//e.encoded[1] = ((ui64)c2 << 32 | e.encoded[1]) / 10;
		$e->encoded[1] = uint32(gmp_div(gmp_or(uint64(gmp_shiftl($c2,32)), $e->encoded[1]),10));
		//unsigned c4 = ((ui64)c3 << 32 | e.encoded[0]) % 10;
		$c4 = uint32(gmp_mod(gmp_or(uint64(gmp_shiftl($c3,32)), $e->encoded[0]),10));
		//e.encoded[0] = ((ui64)c3 << 32 | e.encoded[0]) / 10;
		$e->encoded[0] = uint32(gmp_div(gmp_or(uint64(gmp_shiftl($c3,32)), $e->encoded[0]),10));
		$decimal[34 - $i] = uint8($c4);
	}
	//assert(e.encoded[0] == 0 && e.encoded[1] == 0 && e.encoded[2] == 0 && e.encoded[3] == 0);
	//CHARTYPE* q = confirmation_id;
	$q = $confirmation_id;
	$qi = 0;
	// convert $decimal gmp to int
	for($iii=0;$iii<35;$iii++) {
		$decimal[$iii]=gmp_intval($decimal[$iii]);
	}
	for ($i = 0; $i < 7; $i++) {
		if ($i) { //*q++ = '-';
			$q[$qi] = '-';
			$qi++;
		}
			
		$pi = $i*5;
		/*unsigned char* p = decimal + i*5;
		q[0] = p[0] + '0';
		q[1] = p[1] + '0';
		q[2] = p[2] + '0';
		q[3] = p[3] + '0';
		q[4] = p[4] + '0';
		*/
		$zero = ord('0');
		$q[$qi+0] = chr($decimal[$pi+0] + $zero);
		$q[$qi+1] = chr($decimal[$pi+1] + $zero);
		$q[$qi+2] = chr($decimal[$pi+2] + $zero);
		$q[$qi+3] = chr($decimal[$pi+3] + $zero);
		$q[$qi+4] = chr($decimal[$pi+4] + $zero);

		//q[5] = ((p[0]+p[1]*2+p[2]+p[3]*2+p[4]) % 7) + '0';
		$q[$qi+5] = chr((($decimal[$pi+0]+$decimal[$pi+1]*2+$decimal[$pi+2]+$decimal[$pi+3]*2+$decimal[$pi+4]) % 7) + $zero);
		//q += 6;
		$qi = $qi + 6;
	}
	//*q++ = 0;
	$qi++;
	//$q[$qi] = 0;
	// copy $q to $confirmation_id
	$confirmation_id = $q;
	
	return 0;
}



?>
