;////////////////////////////////////////////////
;// file: matrix_row_ops.asm                   //
;// author: Marijn Verschuren                  //
;// purpose: implement basic row operations    //
;////////////////////////////////////////////////
section .text

global matrix_REF


; args:
; rdi	->	array
; rsi	->	n
; rdx	->	m
matrix_REF:
	ret  ; TODO implement c version and inline row ops!