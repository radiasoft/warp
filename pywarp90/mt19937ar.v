mt19937ar

****** Functions:
init_genrand(s:integer) subroutine # initialize mt(0:N-1) with a seed
init_by_array(init_key:integer,key_length:integer) subroutine # initialize by an array with array-length
genrand_int32() integer function # generates a random number on [0,0xffffffff]-interval
genrand_int31() integer function # generates a random number on [0,0x7fffffff]-interval
genrand_real1() real function # generates a random number on [0,1]-real-interval
genrand_real2() real function # generates a random number on [0,1)-real-interval
genrand_real3() real function # generates a random number on (0,1)-real-interval
genrand_res53() real function # generates a random number on [0,1) with 53-bit resolution
mt_initln() subroutine # initialize large number (over 32-bit constant number)

