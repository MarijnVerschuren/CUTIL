;////////////////////////////////////////////////
;// file: matrix_RREF.asm                      //
;// author: Marijn Verschuren                  //
;// purpose: implement RREF operation          //
;////////////////////////////////////////////////
section .text

global matrix_RREF_64
global matrix_RREF_32


;/*!
; * matrix_RREF.asm
; */
;extern "C" void matrix_RREF_64(f64_t* data, u32_t m, u32_t n);
;extern "C" void matrix_RREF_32(f32_t* data, u32_t m, u32_t n);

; ============================================================
; matrix_RREF_64
;
; arguments:
;   rdi -> pointer to matrix data (double*)
;   rsi -> m (number of rows)
;   rdx -> n (number of columns)   [converted to bytes: n*8]
;
; registers:
;   rax  -> j (row iterator / inner index)
;   rcx  -> destination row pointer / offset
;   r8   -> (8*n - 32) AVX256 loop bound
;   r9   -> i (pivot row / column index)
;   r10  -> row pointer / row offset (varies by loop)
;   r11  -> row index for elimination (i2)
;
;   xmm0 -> pivot value
;   xmm1 -> temporary scalar / vector
;   xmm2 -> temporary scalar / vector
;   xmm3 -> elimination factor (scalar)
;
;   ymm1 -> vector accumulator / destination
;   ymm2 -> source row vector
;   ymm3 -> broadcasted scalar factor
;
; ============================================================
matrix_RREF_64:
    shl rdx, 3                  ; rdx = n * 8 (bytes per row)

%ifdef AVX256
    lea r8, [rdx - 32]           ; r8 = 8*n - 32 (last full YMM chunk)
%endif

    xor r9, r9                  ; i = 0

.REF_loop:                      ; for (i = 0; i < m; i++)
    xorps xmm1, xmm1            ; xmm1 = 0.0 (used for zero compare)

    mov rax, r9                 ; j = i
    mov r10, rdx                ; r10 = 8*n
    imul r10, r9                ; r10 = i * 8*n
    add r10, rdi                ; r10 = &matrix[i][0]

.pivot_loop:                    ; search pivot in column i
    movsd xmm0, [r10 + r9 * 8]    ; xmm0 = matrix[j][i]
    ucomisd xmm0, xmm1          ; compare pivot with 0.0
    jnz .pivot_end              ; found non-zero pivot

    add r10, rdx                ; advance to next row
    inc rax                     ; j++
    cmp rax, rsi
    jl .pivot_loop              ; continue if j < m
    jmp .REF_loop_continue      ; no pivot found in this column

.pivot_end:

    cmp r9, rax                 ; j == i ?
    je .no_row_swap             ; skip swap if same row


; --------------------------------------------------------
; row swap: swap row i with row j
; --------------------------------------------------------
    xor rax, rax
    mov rcx, rdx                ; rcx = 8*n
    imul rcx, r9                ; rcx = i * 8*n
    add rcx, rdi                ; rcx = &matrix[i][0]

%ifdef AVX256
.SIMD_swap_loop:
    cmp rax, r8
    jnl .SISD_swap_loop_cmp

    vmovupd ymm1, [rcx + rax]
    vmovupd ymm2, [r10 + rax]
    vmovupd [rcx + rax], ymm2
    vmovupd [r10 + rax], ymm1

    add rax, 32
    jmp .SIMD_swap_loop
%endif

.SISD_swap_loop:
    movsd xmm1, [rcx + rax]
    movsd xmm2, [r10 + rax]
    movsd [rcx + rax], xmm2
    movsd [r10 + rax], xmm1
    add rax, 8

.SISD_swap_loop_cmp:
    cmp rax, rdx
    jl .SISD_swap_loop

.no_row_swap:


; --------------------------------------------------------
; row scale: scale pivot row to 1
; --------------------------------------------------------
	xor rax, rax

	movsd xmm1, qword [rel f64_one]	; xmm1 = 1.0
	divsd xmm1, xmm0				; xmm1 = 1.0 / xmm0

%ifdef AVX256
	vbroadcastsd ymm0, xmm1	; broadcast scalar xmm1 = x -> ymm0 = [x,x,x,x]
.SIMD_scale_loop:
	cmp rax, r8			; i >= 8*m - 32
	jnl .SISD_scale_loop_cmp		; end SIMD loop if less then 256b is left

	vmovupd	ymm1, [r10 + rax]
	vmulpd	ymm1, ymm1, ymm0
	vmovupd	[r10 + rax], ymm1

	add rax, 32
	jmp .SIMD_scale_loop
%endif

.SISD_scale_loop:
    movsd xmm1, [r10 + rax]
    mulsd xmm1, xmm0    ; xmm1 = xmm1 * xmm0
    movsd [r10 + rax], xmm1
	add rax, 8

.SISD_scale_loop_cmp:
	cmp rax, rdx
	jl .SISD_scale_loop


; --------------------------------------------------------
; row add: add scaled version of row i to j
; --------------------------------------------------------
    xor r11, r11

	; check r11 != i for i = 0 (done later too but pre loop check is added to prevent cache misses)
    cmp r11, r9
    je .REF_row_add_skip
.REF_row_add_loop:              ; eliminate below pivot

    xor rax, rax
    mov rcx, rdx
    imul rcx, r11
    add rcx, rdi                ; rcx = &matrix[i2][0]

    movsd xmm3, [rcx + r9 * 8]    ; load matrix[i2][i]
    xorpd xmm3, [rel f64_neg] ; factor = -factor

%ifdef AVX256
    vbroadcastsd ymm3, xmm3     ; broadcast factor

.SIMD_add_loop:
    cmp rax, r8
    jnl .SISD_add_loop_cmp

    vmovupd ymm1, [rcx + rax]
    vmovupd ymm2, [r10 + rax]
    vfmadd231pd ymm1, ymm2, ymm3 ; row[i2] += row[i] * factor
    vmovupd [rcx + rax], ymm1

    add rax, 32
    jmp .SIMD_add_loop
%endif

.SISD_add_loop:
    movsd xmm1, [rcx + rax]
    movsd xmm2, [r10 + rax]
    mulsd xmm2, xmm3
    addsd xmm1, xmm2
    movsd [rcx + rax], xmm1
    add rax, 8

.SISD_add_loop_cmp:
    cmp rax, rdx
    jl .SISD_add_loop

.REF_row_add_skip:
    inc r11
    cmp r11, r9
    je .REF_row_add_skip
    cmp r11, rsi
    jl .REF_row_add_loop

.REF_loop_continue:
    inc r9
    cmp r9, rsi
    jl .REF_loop
    ret




; ============================================================
; matrix_RREF_32
;
; arguments:
;   rdi -> pointer to matrix data (double*)
;   rsi -> m (number of rows)
;   rdx -> n (number of columns)   [converted to bytes: n*4]
;
; registers:
;   rax  -> j (row iterator / inner index)
;   rcx  -> destination row pointer / offset
;   r8   -> (4*n - 32) AVX256 loop bound
;   r9   -> i (pivot row / column index)
;   r10  -> row pointer / row offset (varies by loop)
;   r11  -> row index for elimination (i2)
;
;   xmm0 -> pivot value
;   xmm1 -> temporary scalar / vector
;   xmm2 -> temporary scalar / vector
;   xmm3 -> elimination factor (scalar)
;
;   ymm1 -> vector accumulator / destination
;   ymm2 -> source row vector
;   ymm3 -> broadcasted scalar factor
;
; ============================================================
matrix_RREF_32:
    shl rdx, 2                  ; rdx = n * 4 (bytes per row)

%ifdef AVX256
    lea r8, [rdx - 32]           ; r8 = 4*n - 32 (last full YMM chunk)
%endif

    xor r9, r9                  ; i = 0

.REF_loop:                      ; for (i = 0; i < m; i++)
    xorps xmm1, xmm1            ; xmm1 = 0.0 (used for zero compare)

    mov rax, r9                 ; j = i
    mov r10, rdx                ; r10 = 4*n
    imul r10, r9                ; r10 = i * 4*n
    add r10, rdi                ; r10 = &matrix[i][0]

.pivot_loop:                    ; search pivot in column i
    movss xmm0, [r10 + r9 * 4]  ; xmm0 = matrix[j][i]
    ucomiss xmm0, xmm1          ; compare pivot with 0.0
    jnz .pivot_end              ; found non-zero pivot

    add r10, rdx                ; advance to next row
    inc rax                     ; j++
    cmp rax, rsi
    jl .pivot_loop              ; continue if j < m
    jmp .REF_loop_continue      ; no pivot found in this column

.pivot_end:

    cmp r9, rax                 ; j == i ?
    je .no_row_swap             ; skip swap if same row


; --------------------------------------------------------
; row swap: swap row i with row j
; --------------------------------------------------------
    xor rax, rax
    mov rcx, rdx                ; rcx = 8*n
    imul rcx, r9                ; rcx = i * 8*n
    add rcx, rdi                ; rcx = &matrix[i][0]

%ifdef AVX256
.SIMD_swap_loop:
    cmp rax, r8
    jnl .SISD_swap_loop_cmp

    vmovups ymm1, [rcx + rax]
    vmovups ymm2, [r10 + rax]
    vmovups [rcx + rax], ymm2
    vmovups [r10 + rax], ymm1

    add rax, 32
    jmp .SIMD_swap_loop
%endif

.SISD_swap_loop:
    movss xmm1, [rcx + rax]
    movss xmm2, [r10 + rax]
    movss [rcx + rax], xmm2
    movss [r10 + rax], xmm1
    add rax, 4

.SISD_swap_loop_cmp:
    cmp rax, rdx
    jl .SISD_swap_loop

.no_row_swap:


; --------------------------------------------------------
; row scale: scale pivot row to 1
; --------------------------------------------------------
	mov     eax, 0x3F800000			; eax = 1.0
	movd	xmm1, eax  				; xmm1 = 1.0
	divss	xmm1, xmm0				; xmm1 = 1.0 / xmm0

	xor rax, rax
%ifdef AVX256
	vbroadcastss ymm0, xmm1	; broadcast scalar xmm1 = x -> ymm0 = [x,x,x,x]
.SIMD_scale_loop:
	cmp rax, r8			; i >= 8*m - 32
	jnl .SISD_scale_loop_cmp		; end SIMD loop if less then 256b is left

	vmovups	ymm1, [r10 + rax]
	vmulps	ymm1, ymm1, ymm0
	vmovups	[r10 + rax], ymm1

	add rax, 32
	jmp .SIMD_scale_loop
%endif

.SISD_scale_loop:
    movss xmm1, [r10 + rax]
    mulss xmm1, xmm0    ; xmm1 = xmm1 * xmm0
    movss [r10 + rax], xmm1
	add rax, 4

.SISD_scale_loop_cmp:
	cmp rax, rdx
	jl .SISD_scale_loop


; --------------------------------------------------------
; row add: add scaled version of row i to j
; --------------------------------------------------------
    xor r11, r11

	; check r11 != i for i = 0 (done later too but pre loop check is added to prevent cache misses)
    cmp r11, r9
    je .REF_row_add_skip
.REF_row_add_loop:              ; eliminate below pivot

    xor rax, rax
    mov rcx, rdx
    imul rcx, r11
    add rcx, rdi                ; rcx = &matrix[i2][0]

    movss xmm3, [rcx + r9 * 4]    ; load matrix[i2][i]
    xorps xmm3, [rel f32_neg]	; factor = -factor

%ifdef AVX256
    vbroadcastss ymm3, xmm3     ; broadcast factor

.SIMD_add_loop:
    cmp rax, r8
    jnl .SISD_add_loop_cmp

    vmovups ymm1, [rcx + rax]
    vmovups ymm2, [r10 + rax]
    vfmadd231ps ymm1, ymm2, ymm3 ; row[i2] += row[i] * factor
    vmovups [rcx + rax], ymm1

    add rax, 32
    jmp .SIMD_add_loop
%endif

.SISD_add_loop:
    movss xmm1, [rcx + rax]
    movss xmm2, [r10 + rax]
    mulss xmm2, xmm3
    addss xmm1, xmm2
    movss [rcx + rax], xmm1
    add rax, 4

.SISD_add_loop_cmp:
    cmp rax, rdx
    jl .SISD_add_loop

.REF_row_add_skip:
    inc r11
    cmp r11, r9
    je .REF_row_add_skip
    cmp r11, rsi
    jl .REF_row_add_loop

.REF_loop_continue:
    inc r9
    cmp r9, rsi
    jl .REF_loop
    ret


; TODO: use immediates and regs!!!
section .rodata
align 16
f64_neg:	dq 0x8000000000000000	; -0.0
f64_one:	dq 0x3FF0000000000000   ;  1.0
f32_neg:	dd 0x80000000			; -0.0