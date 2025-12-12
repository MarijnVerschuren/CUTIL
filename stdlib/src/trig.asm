[bits 64]
section .text
	global sin
	global cos
	global tan
	global csc
	global sec
	global cot


sin:                                    ; double precision sin(x)
	sub rsp, 8                              ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	fld qword [rsp]                         ; load double precision float from the stack to the FPU register
	fsin                                    ; call the sin function from the x87 FPU
	fstp qword [rsp]                        ; load double precision float back to the stack from the FPU register
	movsd xmm0, [rsp]                       ; move sin(x) to xmm0
	add rsp, 8                              ; clean up the stack
	ret

cos:                                    ; double precision cos(x)
	sub rsp, 8                              ; free up some space on the stack
	movsd [rsp], xmm0                       ; store x on the stack
	fld qword [rsp]                         ; load double precision float from the stack to the FPU register
	fcos                                    ; call the cos function from the x87 FPU
	fstp qword [rsp]                        ; load double precision float back to the stack from the FPU register
	movsd xmm0, [rsp]                       ; move cos(x) to xmm0
	add rsp, 8                              ; clean up the stack
	ret

tan:                                    ; double precision tan(x)
	movsd xmm2, xmm0                        ; save x to xmm2
	call cos                                ; calculate cos(x)
	movsd xmm1, xmm0                        ; mov cos(x) to xmm1
	movsd xmm0, xmm2                        ; restore x
	call sin                                ; calculate sin(x)
	divsd xmm0, xmm1                        ; calculating sin(x) / cos(x) => tan(x)
	ret

csc:                                    ; double precision csc(x)
	call sin                                ; calculate cos(x)
	movsd xmm1, xmm0                        ; mov cos(x) to xmm1
	mov rax, 0x3ff0000000000000             ; moving 1.0 into rax because xmm registers cant be set using a literal
	movq xmm0, rax                          ; set xmm0 to 1.0
	divsd xmm0, xmm1                        ; calculating 1 / sin(x) => csc(x)
	ret

sec:                                    ; double precision csc(x)
	call cos                                ; calculate cos(x)
	movsd xmm1, xmm0                        ; mov cos(x) to xmm1
	mov rax, 0x3ff0000000000000             ; moving 1.0 into rax because xmm registers cant be set using a literal
	movq xmm0, rax                          ; set xmm0 to 1.0
	divsd xmm0, xmm1                        ; calculating 1 / sin(x) => csc(x)
	ret

cot:                                    ; double precision tan(x)
	movsd xmm2, xmm0                        ; save x to xmm2
	call sin                                ; calculate sin(x)
	movsd xmm1, xmm0                        ; mov cos(x) to xmm1
	movsd xmm0, xmm2                        ; restore x
	call cos                                ; calculate cos(x)
	divsd xmm0, xmm1                        ; calculating sin(x) / cos(x) => tan(x)
	ret
