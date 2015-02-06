module moleculardynamics
    implicit none
    private
    public sum_forces

contains

    function fill_init_mom(n,momenta) result (momenta)
        mean,standdev=0,5
        momenta(:,0)=np.random.normal(mean,standdev,n)
        momenta(:,1)=np.random.normal(mean,standdev,n)
        momenta(:,2)=np.random.normal(mean,standdev,n)

    function fill_init_pos(m, positions) result(positions)
        positions(0)=(0.0,0.0,0.0)
        positions(1)=(0.5,0.5,0.0)
        positions(2)=(0.5,0.0,0.5)
        positions(3)=(0.0,0.5,0.5)
        integer :: counter=0
        DO i (0,m)
            DO j (0,m)
                DO k (0,m)
                    positions(counter:counter+4)=positions(0:4)+(i,j,k)
                    counter = counter + 4
                end DO
            end DO
        end DO
    end function
    
    function def calc_force(particle1,particle2,m,sigma,epsilon) result([Fx,Fy,Fz])
        integer :: deltax= particle2[0]-particle1[0]
        integer :: delx2= (m-deltax)
        integer :: deltay=(particle2[1]-particle1[1])
        integer :: dely2= (m-deltay)
        integer :: deltaz=(particle2[2]-particle1[2])
        integer :: delz2= (m-deltaz)
        if abs(deltax)>abs(delx2) then
            deltax=delx2
        if abs(deltay)>abs(dely2) then
            deltay=dely2
        if abs(deltaz)>abs(delz2) then
            deltaz=delz2
        r=(deltax**2+deltay**2+deltaz**2)**(1./2)
        F=4*epsilon*((12*sigma**12)/r**13 - (6*sigma**6)/r**7)
        integer :: Fx=F*deltax/r 
        integer :: Fy=F*deltay/r 
        integer :: Fz=F*deltaz/r
    end function

    function changemom(deltat,momenta,forces) result(momenta)
        momenta += forces*deltat
    end function
    
    function changepos(deltat,momenta,positions,mass) result(positions)
        DO i (0,m)
            positions[i] += momenta[i]*(deltat/mass)
            Do j (0,3)
                if positions[j]>m
                    positions[j] = positions[j]-m
           end DO
        end DO
    end function
    
    function sum_forces(positions,forces,Np,m,sigma,epsilon) result(forces)
	Do i (0,Np)
		DO j (0,Np)
			if i/=j then
				forces(i)+= call calc_force(positions(i),positions(j),m,sigma,epsilon) 
		end DO
	end Do
    end function

end module       
        
    
    
 
    
    
       	
        