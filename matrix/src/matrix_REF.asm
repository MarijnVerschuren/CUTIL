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
;
; vars:
; xmm0		->  pivot [0:64]
; xmm1		->  scalar [0:64]
; (xyz)mm2	->	pivot row-let
; (xyz)mm3	->  action row-let
; (xyz)mm4	->	broadcast scalar x4
; rax		->	i		->	for (i = 0; i < n && i < m; i++) {...}
; rcx		->  pivot_col
; r8		->  i2
; r9		->  (sete cmd) and (j)
; r10		->  &array[i * m] -> pivoting_row
; r11		->  8 * m
; ebx		->  &array[i2 * m] -> acting_row
; r12		->  8 * m


matrix_REF:
%ifdef AVX256
    push    rbp                    ; standard prologue
    mov     rbp, rsp
    sub     rsp, 64                ; allocate stack frame (aligned) and scratch

    ; Save callee-saved registers we will use
    mov     [rbp - 8], rbx         ; save rbx
    mov     [rbp -16], r12         ; save r12

    ; Setup locals / easier names
    mov     r10, rdi               ; r10 <- base pointer (A)            ; copy base ptr
    xor     rax, rax               ; rax = i = 0                        ; outer loop index

    mov     r11, rdx               ; r11 = m
    shl     r11, 3                 ; r11 = m * 8 (bytes per row stride) ; compute byte stride

	.outer_loop:
    cmp     rax, rsi               ; if i >= n break
    jge     .done

    ; find pivot_col >= 0 such that A[i*m + pivot_col] != 0
    mov     rcx, rax               ; rcx = pivot_col candidate (start at row i)
	.find_pivot:
    cmp     rcx, rdx               ; if pivot_col >= m -> no pivot in row
    jge     .no_pivot_found

    movsd   xmm0, [r10 + rcx*8]    ; xmm0 = A[ rcx ] (scalar double)    ; load candidate
    ptest   xmm0, xmm0             ; set ZF=1 if xmm0 == 0
    sete    r9b                    ; r9 = 1 if zero, 0 if non-zero      ; note sete sets low byte
    test    r9b, r9b               ; set flags for conditional jump
    jz      .got_pivot             ; if non-zero (ZF==0) we have pivot
    inc     rcx                    ; else pivot candidate ++
    jmp     .find_pivot

	.got_pivot:
    ; At this point:
    ;   rax = i (current row index)
    ;   rcx = pivot_col (index of pivot in this row)
    ;   xmm0 = A[pivot_col] (non-zero pivot scalar)
    ; Prepare to loop over subsequent rows i2 = i+1 .. n-1
    lea     r8, [rax + 1]          ; r8 = i2 start (i+1)
    lea     rbx, [r10 + r11]       ; rbx = pointer to next row (A + stride) ; rbx used as per original

	.row_loop:
    cmp     r8, rsi                ; if i2 >= n done with rows
    jge     .after_rows

    ; Load scalar x = A[i2 * m + pivot_col] ; using rbx as base + rcx*8
    movsd   xmm1, [rbx + rcx*8]    ; xmm1 = scalar from row i2 at column pivot_col
    vdivsd  xmm1, xmm1, xmm0       ; xmm1 = xmm1 / pivot (compute multiplier)

    vbroadcastsd ymm4, xmm1        ; ymm4 = [x,x,x,x] broadcast multiplier

    ; Now apply to columns j = pivot_col .. m-1:
    mov     r12, rcx               ; r12 = j (start)
	.apply_loop:
    cmp     r12, rdx               ; if j >= m break
    jae     .apply_done

    ; process in 4-wide chunks while at least 4 remain
    mov     r9, rdx
    sub     r9, r12                ; r9 = remaining columns
    cmp     r9, 4
    jb      .apply_tail            ; if less than 4 -> handle tail below

    vmovupd ymm2, [r10 + r12*8]    ; ymm2 = A[j..j+3] from base row (unaligned ok)
    vmovupd ymm3, [rbx + r12*8]    ; ymm3 = B[j..j+3] from current row (unaligned ok)
    vmulpd  ymm1, ymm4, ymm2       ; ymm1 = ymm4 * ymm2  (x * A_chunk)
    vsubpd  ymm3, ymm3, ymm1       ; ymm3 = ymm3 - ymm1   (B_chunk -= x*A_chunk)
    vmovupd [rbx + r12*8], ymm3    ; store back to B

    add     r12, 4                 ; j += 4
    jmp     .apply_loop

	.apply_tail:
    ; Handle tail (0..3 remaining columns) with scalar ops to avoid OOB
    cmp     r12, rdx               ; if j >= m nothing to do
    jae     .apply_done

    ; iterate remaining columns one by one
	.tail_scalar_loop:
    movsd   xmm2, [r10 + r12*8]    ; xmm2 = A[j]
    movsd   xmm3, [rbx + r12*8]    ; xmm3 = B[j]
    vmulsd  xmm2, xmm2, xmm1       ; xmm2 = x * A[j]   (note: xmm1 still has scalar multiplier)
    vsubsd  xmm3, xmm3, xmm2       ; xmm3 = B[j] - xmm2
    movsd   [rbx + r12*8], xmm3    ; store scalar back to B[j]
    inc     r12
    cmp     r12, rdx
    jl      .tail_scalar_loop

	.apply_done:
    add     rbx, r11               ; move rbx to next row (rbx += row_stride)
    inc     r8                      ; i2++
    jmp     .row_loop

	.after_rows:
    add     r10, r11               ; advance base ptr to next row (A += stride)
    inc     rax                    ; i++
    jmp     .outer_loop

	.no_pivot_found:
    ; No pivot found in this row (rcx >= m). treat as finished for this row:
    add     r10, r11               ; advance base pointer
    inc     rax
    jmp     .outer_loop

	.done:
    ; restore callee-saved registers and stack
    mov     r12, [rbp -16]         ; restore r12
    mov     rbx, [rbp - 8]         ; restore rbx
    mov     rsp, rbp
    pop     rbp
    ret


%else
%endif




; YMM are volatile
;sub rsp, 80								; free up some space on the stack
;movsd [rsp],		xmm0				; store xmm0[0:64] on the stack		0 - 7
;movsd [rsp + 8],	xmm1				; store xmm1[0:64] on the stack		8 - 15
;vmovupd [rsp + 16],	ymm2				; store ymm2[0:256] on the stack	16 - 47
;vmovupd [rsp + 48],	ymm3				; store ymm3[0:256] on the stack	48 - 79
;
;
;
;vmovupd [rsp + 48],	ymm3				; store ymm3[0:256] on the stack
;vmovupd [rsp + 16],	ymm2				; store ymm2[0:256] on the stack
;movsd [rsp + 8],	xmm1				; restore xmm1[0:64] from the stack
;movsd [rsp],		xmm0				; restore xmm0[0:64] from the stack
;add rsp, 80								; clean up the stack
;ret