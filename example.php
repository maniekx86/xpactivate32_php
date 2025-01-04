<?php

/*
Example use of this code:
We have example Windows XP installation ID: 070846-984606-199354-110900-862045-444713-121295-109842-549701\n");
The result should be: 204455-473963-727110-378854-169136-901266-745572 (result by xp_activate32.exe)

$out = ""; 
$input = "070846-984606-199354-110900-862045-444713-121295-109842-549701";
$result = generate($input, $out); // Result code will be stored to $result and generated activation id to $out

$out -> 204455-473963-727110-378854-169136-901266-745572

Error codes: 
0 is success
ERR_TOO_SHORT = 1
ERR_TOO_LARGE = 2
ERR_INVALID_CHARACTER = 3
ERR_INVALID_CHECK_DIGIT = 4
ERR_UNKNOWN_VERSION = 5
ERR_UNLUCKY = 6;

*/

include_once 'xpgen.php';

$in = $_GET["key"];

$out = "";
$result = generate($in, $out);

if($result==0) {
    echo("OK<br>result: $out");
    die();
} else {
    echo("Failed<br>error code: $result");
}
die();
?>
