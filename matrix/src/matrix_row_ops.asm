;////////////////////////////////////////////////
;// file: matrix_row_ops.asm                   //
;// author: Marijn Verschuren                  //
;// purpose: implement basic row operations    //
;////////////////////////////////////////////////
section .text

global matrix_swap_rows
global matrix_add_rows


; args:
;	rdi		->	data
;	rsi		->	m		-> 8 * m
;	rdx		->	d(est)
;	rcx		->	s(rc)
;
; vars:
;	rax		->	i
;	ymm0	->	dst row-let
;	ymm1	->	src row-let
;	r8		->	8*m - 32 (for AVX256)
matrix_swap_rows:
	xor rax, rax
	shl rsi, 3			; m *= 8
	imul rdx, rsi		; d *= 8*m (d row offset)
	imul rcx, rsi		; s *= 8*m (s row offset)
	add rdx, rdi		; d*m*8 + data (d row ptr)
	add rcx, rdi		; s*m*8 + data (s row ptr)
%ifdef AVX256
	lea r8, [rsi - 32]	; r8 = 8*m - 32
	.SIMD_loop:
	cmp rax, r8			; i >= 8*m - 32
	jnl .SISD_loop		; end SIMD loop if less then 256b is left
	vmovupd ymm0, [rdx + rax]
	vmovupd ymm1, [rcx + rax]
	vmovupd [rdx + rax], ymm1
	vmovupd [rcx + rax], ymm0
	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:
    movsd   xmm0, [rdx + rax]
    movsd   xmm1, [rcx + rax]
    movsd   [rdx + rax], xmm1
    movsd   [rcx + rax], xmm0
	add rax, 8
	cmp rax, rsi
	jl .SISD_loop
	ret



; args:
;	rdi		->	data
;	rsi		->	m		-> 8 * m
;	rdx		->	d(est)
;	rcx		->	s(rc)
;	xmm0	->	scalar	-> [scalar] * 4
;
; vars:
;	rax		->	i
;	ymm1	->	dst row-let
;	ymm2	->	src row-let
;	r8		->	8*m - 32 (for AVX256)
matrix_add_rows:
	xor rax, rax
	shl rsi, 3				; m *= 8
	imul rdx, rsi			; d *= 8*m (d row offset)
	imul rcx, rsi			; s *= 8*m (s row offset)
	add rdx, rdi			; d*m*8 + data (d row ptr)
	add rcx, rdi			; s*m*8 + data (s row ptr)
%ifdef AVX256
	vbroadcastsd ymm0, xmm0	; broadcast scalar xmm0 = x -> ymm0 = [x,x,x,x]
	lea r8, [rsi - 32]		; r8 = 8*m - 32
	.SIMD_loop:
	cmp rax, r8				; i >= 8*m - 32
	jnl .SISD_loop			; end SIMD loop if less then 256b is left

	vmovupd ymm1, [rdx + rax]
	vmovupd ymm2, [rcx + rax]
	vfmadd231pd ymm1, ymm2, ymm0	; ymm1 = ymm1 + ymm2 * ymm0
	vmovupd [rdx + rax], ymm1

	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:

    movsd xmm1, [rdx + rax]
    movsd xmm2, [rcx + rax]
    mulsd xmm2, xmm0    ; xmm2 = xmm2 * xmm0
	addsd xmm1, xmm2    ; xmm1 += xmm2
    movsd [rdx + rax], xmm1

	add rax, 8
	cmp rax, rsi
	jl .SISD_loop
	ret


