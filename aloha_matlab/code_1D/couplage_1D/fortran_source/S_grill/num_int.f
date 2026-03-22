			
c polynom interpolation
			SUBROUTINE polint(xa,ya,n,x,y,dy)
			INTEGER n,NMAX
			REAL dy,x,y,xa(n),ya(n)
c			Largest anticipated value of n.			
			PARAMETER (NMAX=10) 
c			Given arrays xa and ya, each of length n, and given a value x, this routine returns a
c			value y, and an error estimate dy. If P(x) is the polynomial of degree N ? 1 such that
c			P(xai) = yai, i = 1, . . . , n, then the returned value y = P(x).
			INTEGER i,m,ns
			REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
			ns=1
			dif=abs(x-xa(1))
c			Here we find the index ns of the closest table entry,
			do i=1,n 
				dift=abs(x-xa(i))
				if (dift.lt.dif) then
					ns=i
					dif=dift
				endif
c			and initialize the tableau of c?s and d?s.
				c(i)=ya(i) 
				d(i)=ya(i)
			enddo
c			This is the initial approximation to y.
			y=ya(ns) 
			ns=ns-1
c			For each column of the tableau,
			do m=1,n-1 
c			we loop over the current c?s and d?s and update them.
				do i=1,n-m 
					ho=xa(i)-x
					hp=xa(i+m)-x
					w=c(i+1)-d(i)
					den=ho-hp
					if(den.eq.0.) write (*,*) 'failure in polint'
c			This error can occur only if two input xa?s are (to within roundoff) identical.
					den=w/den
c			Here the c?s and d?s are updated.
					d(i)=hp*den 
					c(i)=ho*den
				enddo
				if (2*ns.lt.n-m) then 
c			After each column in the tableau is completed, we decide
c			which correction, c or d, we want to add to our accumulating
c			value of y, i.e., which path to take through
c			the tableau?forking up or down. We do this in such a
c			way as to take the most ?straight line? route through the
c			tableau to its apex, updating ns accordingly to keep track
c			of where we are. This route keeps the partial approximations
c			centered (insofar as possible) on the target x. T he
c			last dy added is thus the error indication.
					dy=c(ns+1)
				else
					dy=d(ns)
					ns=ns-1
				endif
				y=y+dy
			enddo
			return
			END


c integration par la methode de Romberg		
			SUBROUTINE qromb(func,a,b,ss)
			INTEGER JMAX,JMAXP,K,KM
			REAL a,b,func,ss,EPS
			EXTERNAL func
			PARAMETER (EPS=1.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
c			USES polint,trapzd
c			Returns as ss the integral of the function func from a to b. Integration is performed by
c			Romberg's method of order 2K, where, e.g., K=2 is Simpson?s rule.
c			Parameters: EPS is the fractional accuracy desired, as determined by the extrapolation
c			error estimate; JMAX limits the total number of steps; K is the number of points used in
c			the extrapolation.
			INTEGER j
c			These store the successive trapezoidal approximations
			REAL dss,h(JMAXP),s(JMAXP) 
c			and their relative stepsizes.
			h(1)=1. 
			do j=1,JMAX
				call trapzd(func,a,b,s(j),j)
				if (j.ge.K) then
					call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
					if (abs(dss).le.EPS*abs(ss)) return
				endif
				s(j+1)=s(j)
			
c			This is a key step: The factor is 0.25 even though
c			the stepsize is decreased by only 0.5. This makes
c			the extrapolation a polynomial in h2 as allowed
c			by equation (4.2.1), not just a polynomial in h.

				h(j+1)=0.25*h(j) 
			enddo
c			write (*,*) 'too many steps in qromb'
			END		

c subroutine d'integration classique (trapeze)
			SUBROUTINE trapzd(func,a,b,s,n)
			INTEGER n
			REAL a,b,s,func
			EXTERNAL func
c This routine computes the nth stage of refinement of an extended trapezoidal rule. func is
c input as the name of the function to be integrated between limits a and b, also input. When
c called with n=1, the routine returns as s the crudest estimate of integral from a to 0 of f(x)dx. 
c Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2^(n-2)
c additional interior points. s should not be modified between sequential calls.			
			INTEGER it,j
			REAL del,sum,tnm,x
			if (n.eq.1) then
				s=0.5*(b-a)*(func(a)+func(b))
			else
				it=2**(n-2)
				tnm=it
c				This is the spacing of the points to be added.				
				del=(b-a)/tnm 
				x=a+0.5*del
				sum=0.
				do j=1,it
					sum=sum+func(x)
					x=x+del
				enddo
c				This replaces s by its refined value.
				s=0.5*(s+(b-a)*sum/tnm) 
			endif
			return
			END

c *********************************************************
			SUBROUTINE qtrap(func,a,b,s)
			INTEGER JMAX
			REAL a,b,func,s,EPS
			EXTERNAL func
c MODIF JH : JMAX=12 au lieu de 20 (pour acc�l�rer les calculs au d�triment de la convergence)			
			PARAMETER (EPS=1.e-6, JMAX=20)
c USES trapzd
c			Returns as s the integral of the function func from a to b. The parameters EPS can be set
c			to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum
c			allowed number of steps. Integration is performed by the trapezoidal rule.
			INTEGER j
			REAL olds
c 		Any number that is unlikely to be the average of the function
c			at its endpoints will do here.
			olds=-1.e30 
			do j=1,JMAX 
				call trapzd(func,a,b,s,j)
c			Avoid spurious early convergence.
				if (j.gt.5) then 
					if (abs(s-olds).lt.EPS*abs(olds).or.(s.eq.0.and.olds.eq.0.)) return
				endif
				olds=s
			enddo
c			pause 'too many steps in qtrap'
c			write (*,*) '[S_grill_version6] (WARNING) : too many steps in qtrap'
			END
