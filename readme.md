# xpactivate32_php
PHP port of xp_activate32.exe

This is a full PHP port of xp_activate32.exe. Ported line by line. But why? - because why not! And to be honest, I don't understand any of the calculations, code was made in a way to replicate the original code workings in 1:1.

Since it's a "1:1 port", before submitting a issue with an installation id that doesn't work, check if it works in the original xp\_activate32.exe. If it doesn't work in this port, but works with xp\_activate32.exe - then go ahead. 

Source used: https://archive.org/details/xp_activate32_src

# Usage
xpgen.php is the main file.
```php
$out = ""; 
$input = "070846-984606-199354-110900-862045-444713-121295-109842-549701";
$result = generate($input, $out); // Result code will be stored to $result and generated activation id to $out

// $out -> 204455-473963-727110-378854-169136-901266-745572
```

## Error codes: 
```
SUCCESS = 0
ERR_TOO_SHORT = 1
ERR_TOO_LARGE = 2
ERR_INVALID_CHARACTER = 3
ERR_INVALID_CHECK_DIGIT = 4
ERR_UNKNOWN_VERSION = 5
ERR_UNLUCKY = 6;
```

# Live demo
A working live demo of this port with a simple frontend is available [here].(https://maniek86.xyz/dev/xpgen/xpgen.html). 
