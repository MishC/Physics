program divisors
!	This program finds the divisors of an integer input by the user.
!	The divisors are printed to a file.
integer n, i, j, k, d(10)
open (unit = 1, file = "divisorses.dat")
print *, "Enter a positive integer :"
read *, n
write (1,*) "Here are the divisors of ", n, " :"
k = 0
do i = 1, n
if (mod(n,i) .eq. 0) then
k = k + 1
d(k) = i
end if
if (k .eq. 10) then
write (1,5) (d(j), j = 1, 10)
k = 0
end if
end do
write (1,5) (d(j), j = 1, k)
  5	format (10i7)
close (1)
print *, "The divisors are listed in the file 'divisors'. Bye."
end
