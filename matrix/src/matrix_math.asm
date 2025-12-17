;////////////////////////////////////////////////
;// file: matrix_math.asm                      //
;// author: Marijn Verschuren                  //
;// purpose: implement REF operation           //
;////////////////////////////////////////////////
section .text

global matrix_add_64
global matrix_add_32
global matrix_sub_64
global matrix_sub_32
global matrix_scale_64
global matrix_scale_32


; ============================================================
; matrix_add/sub_64/32
;
; arguments:
;   rdi	->	dst
;   rsi	->	a
;   rdx	->	b
;   rcx	->	cnt		-> 8*cnt
;
; vars:
;	rax	->	i
;	r8	->	8*cnt - 32
; ============================================================
matrix_add_64:
	shl rcx, 3		; cnt *= 8
	xor rax, rax	; i = 0
%ifdef AVX256
	lea r8, [rcx - 32]		; r8 = 8*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 8*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovupd ymm0, [rsi + rax]	; load ymm0 = [a0,a1,a2,a3] starting from rsi + i
	vmovupd ymm1, [rdx + rax]	; load ymm1 = [b0,b1,b2,b3] starting from rdx + i
	vaddpd ymm0, ymm0, ymm1		; ymm0 += ymm1
	vmovupd [rdi + rax], ymm0	; store [a0 + b0, a1 + b1, a2 + b2, a3 + b3]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movsd xmm0, [rsi + rax]	; load xmm0 = a at rsi + i
    movsd xmm1, [rdx + rax]	; load xmm1 = b at rdx + i
	addsd xmm0, xmm1    ; xmm0 += xmm1
    movsd [rdi + rax], xmm0

	add rax, 8
	.SISD_loop_cmp
	cmp rax, rcx
	jl .SISD_loop
	ret


matrix_add_32:
	shl rcx, 2		; cnt *= 4
	xor rax, rax	; i = 0
%ifdef AVX256
	lea r8, [rcx - 32]		; r8 = 4*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 4*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovups ymm0, [rsi + rax]	; load ymm0 = [a0,... ,a7] starting from rsi + i
	vmovups ymm1, [rdx + rax]	; load ymm1 = [b0,... ,b7] starting from rdx + i
	vaddps ymm0, ymm0, ymm1		; ymm0 += ymm1
	vmovups [rdi + rax], ymm0	; store [a0 + b0,... , a7 + b7]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movss xmm0, [rsi + rax]	; load xmm0 = a at rsi + i
    movss xmm1, [rdx + rax]	; load xmm1 = b at rdx + i
	addss xmm0, xmm1    ; xmm0 += xmm1
    movss [rdi + rax], xmm0

	add rax, 8
	.SISD_loop_cmp
	cmp rax, rcx
	jl .SISD_loop
	ret


matrix_sub_64:
	shl rcx, 3		; cnt *= 8
	xor rax, rax	; i = 0
%ifdef AVX256
	lea r8, [rcx - 32]		; r8 = 8*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 8*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovupd ymm0, [rsi + rax]	; load ymm0 = [a0,a1,a2,a3] starting from rsi + i
	vmovupd ymm1, [rdx + rax]	; load ymm1 = [b0,b1,b2,b3] starting from rdx + i
	vsubpd ymm0, ymm0, ymm1		; ymm0 -= ymm1
	vmovupd [rdi + rax], ymm0	; store [a0 + b0, a1 + b1, a2 + b2, a3 + b3]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movsd xmm0, [rsi + rax]	; load xmm0 = a at rsi + i
    movsd xmm1, [rdx + rax]	; load xmm1 = b at rdx + i
	subsd xmm0, xmm1    ; xmm0 -= xmm1
    movsd [rdi + rax], xmm0

	add rax, 8
	.SISD_loop_cmp
	cmp rax, rcx
	jl .SISD_loop
	ret


matrix_sub_32:
	shl rcx, 2		; cnt *= 4
	xor rax, rax	; i = 0
%ifdef AVX256
	lea r8, [rcx - 32]		; r8 = 4*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 4*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovups ymm0, [rsi + rax]	; load ymm0 = [a0,... ,a7] starting from rsi + i
	vmovups ymm1, [rdx + rax]	; load ymm1 = [b0,... ,b7] starting from rdx + i
	vsubps ymm0, ymm0, ymm1		; ymm0 -= ymm1
	vmovups [rdi + rax], ymm0	; store [a0 + b0,... , a7 + b7]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movss xmm0, [rsi + rax]	; load xmm0 = a at rsi + i
    movss xmm1, [rdx + rax]	; load xmm1 = b at rdx + i
	subss xmm0, xmm1    ; xmm0 -= xmm1
    movss [rdi + rax], xmm0

	add rax, 8
	.SISD_loop_cmp
	cmp rax, rcx
	jl .SISD_loop
	ret




; ============================================================
; matrix_scale_64/32
;
; arguments:
;   rdi		->	dst
;   rsi		->	src
;   rdx		->	cnt		-> 8*cnt
;   xmm0	->	scalar
;
; vars:
;	rax	->	i
;	r8	->	8*cnt - 32
; ============================================================
matrix_scale_64:
	shl rdx, 3		; cnt *= 8
	xor rax, rax	; i = 0
%ifdef AVX256
	vbroadcastsd ymm0, xmm0	; broadcast scalar xmm0 = x -> ymm0 = [x,x,x,x]
	lea r8, [rdx - 32]		; r8 = 8*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 8*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovupd ymm1, [rsi + rax]	; load ymm1 = [a0,a1,a2,a3] starting from rsi + i
	vmulpd ymm1, ymm1, ymm0
	vmovupd [rdi + rax], ymm1	; store [a0 + b0, a1 + b1, a2 + b2, a3 + b3]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movsd xmm1, [rsi + rax]	; load xmm1 = a at rsi + i
	mulsd xmm1, xmm0	; xmm1 *= xmm0
    movsd [rdi + rax], xmm1

	add rax, 8
	.SISD_loop_cmp
	cmp rax, rdx
	jl .SISD_loop
	ret


matrix_scale_32:
	shl rdx, 2		; cnt *= 4
	xor rax, rax	; i = 0
%ifdef AVX256
	vbroadcastss ymm0, xmm0	; broadcast scalar xmm0 = x -> ymm0 = [x0 ,... ,x7]
	lea r8, [rdx - 32]		; r8 = 4*cnt - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 4*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left

	vmovups ymm1, [rsi + rax]	; load ymm1 = [a0,a1,a2,a3] starting from rsi + i
	vmulps ymm1, ymm1, ymm0
	vmovups [rdi + rax], ymm1	; store [a0 + b0, a1 + b1, a2 + b2, a3 + b3]

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movss xmm1, [rsi + rax]	; load xmm1 = a at rsi + i
	mulss xmm1, xmm0	; xmm1 *= xmm0
    movss [rdi + rax], xmm1

	add rax, 4
	.SISD_loop_cmp
	cmp rax, rdx
	jl .SISD_loop
	ret