bits 64
section .text
	global test


test:
%ifdef __AVX2__
	jmp $
	ret
%else
	ret
%endif