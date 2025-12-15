;////////////////////////////////////////////////
;// file: matrix_row_ops.asm                   //
;// author: Marijn Verschuren                  //
;// purpose: implement basic row operations    //
;////////////////////////////////////////////////
section .text

global matrix_REF


; args:
;	rdi		->	array
;	rsi		->	n
;	rdx		->	m	-> 8m
;
; vars:
;	REF:
;		r9		->	i
;		r10		->	j * 8m
;		xmm0	->	pivot
;	row_ops
;		rax		->	i
;		ymm1	->	row_op_reg
;		ymm2	->	row_op_reg
;		r11		->	row_add_dst
;		r8		->	8*m - 32 (for AVX256)
matrix_REF:
	shl rdx, 3	; m *= 8
%ifdef AVX256
	lea r8, [rdx - 32]	; r8 = 8*m - 32
%endif
	xor r9, r9
	.REF_loop:;	for (uint32_t j, i = 0; i < n; i++)
	;xorps xmm1, xmm1	; xmm1 = 0 TODO: other reg?
	mov rax, r9			; j = i			(indexing variable)
	mov r10, rdx		; j = 8m
	imul r10, r9		; j *= i	(offset ptr to element)
	add r10, rdi		; j += data
	.pivot_loop:; for (j = i; j < n; j++)
    movsd xmm0, [r10 + r9 * 8]	; pivot = data[j * m + i]
	ptest xmm0, xmm0	; pivot == 0.0?
	jnz .pivot_end		; found pivot, i.e pivot != 0.0
	add r10, rdx		; 8mj += 8m (j++)
	inc rax				; j++
	cmp rax, rsi		; j < n
	jl .pivot_loop		; loop if i < n
	jmp .REF_loop_continue	; skip column if no pivots found
	.pivot_end			; found pivot
	cmp r9, rax			; j == i?
	je .no_row_swap		; skip row swap if j == i

	;/*!
	; * swap rows
	; */
	; args:
	;	rdi		->	array
	;	rdx		->	8m
	;	rax		->	i
	; 	r9		->	d(est)	->	REF(i)
	;	r10		->	s(rc)	->	8mj
	;	rcx		->	d(rc) offset -> 8mi
	xor rax, rax
	mov rcx, rdx		; s = 8*m
	imul rcx, r9		; s *= i (s row offset)
	add r10, rdi		; s*m*8 + array (d row ptr)
	add rcx, rdi		; d*m*8 + array (s row ptr)
%ifdef AVX256
	.SIMD_swap_loop:
	cmp rax, r8			; i >= 8*m - 32
	jnl .SISD_swap_loop_cmp		; end SIMD loop if less then 256b is left
	vmovupd ymm1, [rcx + rax]
	vmovupd ymm2, [r10 + rax]
	vmovupd [rcx + rax], ymm2
	vmovupd [r10 + rax], ymm1
	add rax, 32
	jmp .SIMD_swap_loop
%endif
	.SISD_swap_loop:
    movsd   xmm1, [rcx + rax]
    movsd   xmm2, [r10 + rax]
    movsd   [rcx + rax], xmm2
    movsd   [r10 + rax], xmm1
	add rax, 8
    .SISD_swap_loop_cmp:
	cmp rax, rdx
	jl .SISD_swap_loop

	.no_row_swap:

	lea r11, [r9 + 1]
	mov r10, rdx		; s = 8m
	imul r10, r9		; s *= i (s row offset)
	add r10, rdi		; s*m*8 + data (s row ptr)
	.REF_row_add_loop:

	;/*!
	; * add rows
	; */
	xor rax, rax
	mov rcx, rdx		; d = 8m
	imul rcx, r11		; d = j*8*m (d row offset)
	add rcx, rdi		; d*m*8 + data (d row ptr)
	movsd xmm3, [rcx + r9 * 8]
	divsd xmm3, xmm0				; xmm3 = xmm3 / xmm0
	xorpd xmm3, [rel f64_sign_mask]	; xmm3 = -(xmm3 / xmm0)
%ifdef AVX256
	vbroadcastsd ymm3, xmm3	; broadcast scalar xmm0 = x -> ymm0 = [x,x,x,x]
	.SIMD_add_loop:
	cmp rax, r8				; i >= 8*m - 32
	jnl .SISD_add_loop_cmp		; end SIMD loop if less then 256b is left
	vmovupd ymm1, [rcx + rax]
	vmovupd ymm2, [r10 + rax]
	vfmadd231pd ymm1, ymm2, ymm3	; ymm1 = ymm1 + ymm2 * ymm3
	vmovupd [rcx + rax], ymm1
	add rax, 32
	jmp .SIMD_add_loop
%endif
	.SISD_add_loop:
    movsd xmm1, [rcx + rax]
    movsd xmm2, [r10 + rax]
    mulsd xmm2, xmm3    ; xmm2 = xmm2 * xmm3
	addsd xmm1, xmm2    ; xmm1 += xmm2
    movsd [rcx + rax], xmm1
	add rax, 8
    .SISD_add_loop_cmp:
	cmp rax, rdx
	jl .SISD_add_loop

	inc r11
	cmp r11, rsi
	jl .REF_row_add_loop

	.REF_loop_continue:
	inc r9			; i++
	lea r10, [r9 + 1]
	cmp r10, rsi	; (i + 1) < n
	jl .REF_loop	; loop if i < n
	ret


section .rodata
align 16
f64_sign_mask: dq 0x8000000000000000