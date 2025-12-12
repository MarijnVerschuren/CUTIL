[bits 64]
section .text
	global sqrt
	global sqrt32
	global ln
	global log2
	global log10
	global exp
	global pow


sqrt:
	sqrtsd xmm0, xmm0
	ret

ln:
	sub rsp, 16                             ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	mov rax, 0x3fe62e42fefa3bdc             ; move log(2) / log(e) into rax to convert log(x) / log(2) to log(x) / log(e)
	mov [rsp + 8], rax                      ; move log(2) / log(e) onto the stack
	fld qword [rsp + 8]                     ; load double precision float for 1.0 from the stack to the FPU register
	fld qword [rsp]                         ; load double precision float for x from the stack to the FPU register
	fyl2x                                   ; call the log-ab function from the x87 FPU
	fstp qword [rsp]                        ; load double precision float back to the stack from the FPU register
	movsd xmm0, [rsp]                       ; move ln(x) to xmm0
	add rsp, 16                             ; clean up the stack
	ret

log2:
	sub rsp, 8                              ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	fld1                                    ; load double precision float for 1.0 from the stack to the FPU register
	fld qword [rsp]                         ; load double precision float for x from the stack to the FPU register
	fyl2x                                   ; call the log-ab function from the x87 FPU
	fstp qword [rsp]                        ; load double precision float back to the stack from the FPU register
	movsd xmm0, [rsp]                       ; move log2(x) to xmm0
	add rsp, 8                              ; clean up the stack
	ret

log10:
	sub rsp, 16                             ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	mov rax, 0x3fd34413509f7b52             ; move log(2) / log(10) into rax to convert log(x) / log(2) to log(x) / log(10)
	mov [rsp + 8], rax                      ; move log(2) / log(10) onto the stack
	fld qword [rsp + 8]                     ; load double precision float for 1.0 from the stack to the FPU register
	fld qword [rsp]                         ; load double precision float for x from the stack to the FPU register
	fyl2x                                   ; call the log-ab function from the x87 FPU
	fstp qword [rsp]                        ; load double precision float back to the stack from the FPU register
	movsd xmm0, [rsp]                       ; move log10(x) to xmm0
	add rsp, 16                             ; clean up the stack
	ret

exp:
	sub rsp, 8                              ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	fld qword [rsp]                         ; load the current approximation to the FPU register
	fldl2e                                  ; load log2(e) into the FPU register
	fmulp                                   ; multiply x by log2(e)
	f2xm1                                   ; calculate 2 ^ (x * log2(e)) - 1
	fld1                                    ; load 1.0 into the FPU register
	faddp                                   ; add 1.0 to the result of the previous operation
	fstp qword [rsp]                        ; store the result back on the stack
	movsd xmm0, [rsp]                       ; load e^x into xmm0
	add rsp, 8                              ; clean up space on the stack
	ret

pow:
	sub rsp, 8                              ; free up some space on the stack
	call log2                               ; calculate log2(x)
	mulsd xmm0, xmm1                        ; calculate y * log2(x)
	movsd [rsp], xmm0                       ; store y * log2(x) on the stack
	fld qword [rsp]                         ; load the current approximation to the FPU register
	f2xm1                                   ; calculate 2 ^ (y * log2(x)) - 1
	fld1                                    ; load 1.0 into the FPU register
	faddp                                   ; add 1.0 to the result of the previous operation
	fstp qword [rsp]                        ; store the result back on the stack
	movsd xmm0, [rsp]                       ; load x^y into xmm0
	add rsp, 8                              ; clean up space on the stack
	ret
