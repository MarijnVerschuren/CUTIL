;////////////////////////////////////////////////
;// file: matrix_row_ops.asm                   //
;// author: Marijn Verschuren                  //
;// purpose: implement basic row operations    //
;////////////////////////////////////////////////
section .text

global matrix_swap_rows
global matrix_add_rows
global matrix_scale_row


;/*!
; * matrix_row_ops.asm
; */
;extern "C" void matrix_swap_rows(f64_t* data, u32_t m, u32_t d, u32_t s);
;extern "C" void matrix_add_rows(f64_t* data, u32_t m, u32_t d, u32_t s, f64_t scalar);
;extern "C" void matrix_scale_row(f64_t* data, u32_t m, u32_t r, f64_t scalar);


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
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left
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
    .SISD_loop_cmp
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
	jnl .SISD_loop_cmp			; end SIMD loop if less then 256b is left

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
	.SISD_loop_cmp
	cmp rax, rsi
	jl .SISD_loop
	ret


; args:
;	rdi		->	data
;	rsi		->	m		-> 8 * m
;	rdx		->	r(ow)
;	xmm0	->	scalar
;
; vars:
;	rax		->	i
;	ymm1	->	src row-let
;	r8		->	8*m - 32 (for AVX256)
matrix_scale_row:
	xor rax, rax
	shl rsi, 3			; m *= 8
	imul rdx, rsi		; r *= 8*m (row offset)
	add rdx, rdi		; r*m*8 + data (row ptr)
%ifdef AVX256
	vbroadcastsd ymm0, xmm0	; broadcast scalar xmm0 = x -> ymm0 = [x,x,x,x]
	lea r8, [rsi - 32]	; r8 = 8*m - 32
	.SIMD_loop:
	cmp rax, r8			; i >= 8*m - 32
	jnl .SISD_loop_cmp		; end SIMD loop if less then 256b is left
	vmovupd ymm1, [rdx + rax]
	vmulpd ymm1, ymm1, ymm0
	vmovupd [rdx + rax], ymm1
	add rax, 32
	jmp .SIMD_loop
%endif
	.SISD_loop:
    movsd xmm1, [rdx + rax]
    mulsd xmm1, xmm0    ; xmm1 = xmm1 * xmm0
    movsd [rdx + rax], xmm1
	add rax, 8
    .SISD_loop_cmp
	cmp rax, rsi
	jl .SISD_loop
	ret
