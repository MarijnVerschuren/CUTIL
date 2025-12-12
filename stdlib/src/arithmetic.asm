[bits 64]
section .text
	global fabs
	global fabs32


fabs:
	mov rax, 0x7fffffffffffffff             ; move the sign mask into rax
	movq xmm1, rax                          ; move the sign mask into xmm1
	andps xmm0, xmm1                        ; disable the sign bit
	ret
fabs32:
	mov eax, 0x7fffffff                     ; move the sign mask into eax
	movd xmm1, eax                          ; move the sign mask into xmm1
	andps xmm0, xmm1                        ; disable the sign bit
	ret
