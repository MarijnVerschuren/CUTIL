;////////////////////////////////////////////////
;// file: matrix_REF.asm                       //
;// author: Marijn Verschuren                  //
;// purpose: implement REF operation           //
;////////////////////////////////////////////////
section .text

global matrix_REF_64
global matrix_REF_32



; ============================================================
; matrix_REF_64
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
matrix_REF_64:
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
; row add: add scaled version of row i to j
; --------------------------------------------------------
    lea r11, [r9 + 1]           ; i2 = i + 1

.REF_row_add_loop:              ; eliminate below pivot

    xor rax, rax
    mov rcx, rdx
    imul rcx, r11
    add rcx, rdi                ; rcx = &matrix[i2][0]

    movsd xmm3, [rcx + r9 * 8]    ; load matrix[i2][i]
    divsd xmm3, xmm0            ; factor = A[i2][i] / pivot
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

    inc r11
    cmp r11, rsi
    jl .REF_row_add_loop

.REF_loop_continue:
    inc r9
    lea r10, [r9 + 1]
    cmp r10, rsi
    jl .REF_loop
    ret



; ============================================================
; matrix_REF_32
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
matrix_REF_32:
    shl rdx, 2                  ; rdx = n * 4 (bytes per row)

%ifdef AVX256
    lea r8, [rdx - 32]           ; r8 = 4*n - 32 (last full YMM chunk)
%endif

    xor r9, r9                  ; i = 0

.REF_loop:                      ; for (i = 0; i < m; i++)
    xorps xmm1, xmm1            ; xmm1 = 0.0 (used for zero compare)

    mov rax, r9                 ; j = i
    mov r10, rdx                ; r10 = 8*n
    imul r10, r9                ; r10 = i * 8*n
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
    mov rcx, rdx                ; rcx = 4*n
    imul rcx, r9                ; rcx = i * 4*n
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
; row add: add scaled version of row i to j
; --------------------------------------------------------
    lea r11, [r9 + 1]           ; i2 = i + 1

.REF_row_add_loop:              ; eliminate below pivot

    xor rax, rax
    mov rcx, rdx
    imul rcx, r11
    add rcx, rdi                ; rcx = &matrix[i2][0]

    movss xmm3, [rcx + r9 * 4]    ; load matrix[i2][i]
    divss xmm3, xmm0            ; factor = A[i2][i] / pivot
    xorps xmm3, [rel f32_neg] ; factor = -factor

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

    inc r11
    cmp r11, rsi
    jl .REF_row_add_loop

.REF_loop_continue:
    inc r9
    lea r10, [r9 + 1]
    cmp r10, rsi
    jl .REF_loop
    ret



; TODO: use immediates and regs!!!
section .rodata
align 16
f64_neg:	dq 0x8000000000000000	; -0.0
f32_neg:	dd 0x80000000			; -0.0