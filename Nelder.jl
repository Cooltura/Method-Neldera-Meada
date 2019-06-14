
#f(x,y)=(2*x+y*y)*(MathConstants.e)^x
f(x,y)=x*x+y*y
#f(x)=x*x

#P=[[0.0 -0.5],[ 1.0 0.5],[0.5 1.0]]
P=[[0.0 -0.5 0.0],[ 0.5 0.5 1.0],[0.5 1.0 -0.5],[0.5 -1.0 -0.5]]
#P=[7.5 2.0]
length(P)*1.0


function Change(P)
P=[[P[1][1]+0.001 P[1][2]+0.001],[P[2][1]-0.001 P[2][2]-0.001],[P[3][1]+0.001 P[3][2]+0.001]]  
    return P
end

function Center(P,h,n)
    s=zeros(Float64,1,n-1)
    Pc=zeros(Float64,1,n-1)
    s2=0.0
   
         
     for i in 1:n
        if i != h
            for j in 1:n-1
                s[j]+=P[i][j]
             end
            end
        end
    for i in 1:n-1
        Pc[i]=s[i]/n
    end
    return Pc
    
end

function maxx(P,f,n)
    arg_f=zeros(Float64,1,n-1)
    arg_h=zeros(Float64,1,n-1)
    h=1
        for i in 2:n
            for j in 1:n-1
                arg_f[j]=P[i][j] 
                arg_h[j]=P[h][j]
            end
            if f(arg_f...)>f(arg_h...) 
              h=i
            end
        end 
    return h
end

function newmax(P,Pr,f,h,n)
    arg_f=zeros(Float64,1,n-1)
    arg_h=zeros(Float64,1,n-1)
    arg_pr=zeros(Float64,1,n-1)
    if h!=1
      h2=1
    else 
       h2=2
    end
  
         for i in 1:n
            for j in 1:n-1
                arg_f[j]=P[i][j] 
                arg_h[j]=P[h2][j]
                arg_pr[j]=Pr[j]
            end
            if n!=h
              if f(arg_f...)>f(arg_h...) 
               h2=i
              end
            end
          end
            if f(arg_pr...)>f(arg_h...)
             return true
            else 
             return false
            end    
           
end

function min(P,f,n)
    arg_f=zeros(Float64,1,n-1)
    arg_h=zeros(Float64,1,n-1)
    h=1
        for i in 2:n
            for j in 1:n-1
                arg_f[j]=P[i][j]
                arg_h[j]=P[h][j]
            end
           if f(arg_f...)<f(arg_h...) 
            h=i
            end
        end 
    return h
end

function reduction(P,l,n)    
    PL=P[l]
   
        for i in 1:n
            for j in 1:n-1
            P[i][j]=0.5*(P[i][j]+PL[j])
            end
        end 
end


function reflection(Ph,Pc,P,n)
    a=1.0
    Pr=zeros(Float64,1,n-1)
    for i in 1:n-1
    Pr[i]=(1.0+a)*Pc[i]-a*Ph[i]  
    end
    
    return Pr
end

function expansion(Pr,Pc,P,n)
    c=2.0
 
    Pe=zeros(Float64,1,n-1)
    for i in 1:n-1
      Pe[i]=(1.0+c)*Pr[i]-c*Pc[i]
    end
  
    return Pe   
end

function contraction(Ph,Pc,n)
    b=0.5
     Pcon=zeros(Float64,1,n-1)
    for i in 1:n-1
      Pcon[i]=b*Ph[i]+(1.0-b)*Pc[i]
    end
    return Pcon   
end

function Check(P2, P, x,f,n)

    arg_p=zeros(Float64,1,n-1)
    arg_x=zeros(Float64,1,n-1)
     
        for i in 1:n-1
            arg_p[i]=P2[i]
            arg_x[i]=P[x][i]
        end
       if f(arg_p...)<f(arg_x...)
            return true
       else
            return false
       end
end

function count_e(P,n)

    wynik= P[1][1]
    for i in 2:n
        wynik=wynik-P[i][1]
        end
    wynik=wynik*wynik
    return wynik
    
end

function create_P(n)
  if n==1
       P=[[1.0],[2.0]]
  else
   P = [rand(n) for i in 1:n+1]
   # P=vec(P)
   end
        return P
end


function m_Nelder(f,n)
    n=n+1
step=0
t=0.00000001
e=1.0
P=create_P(n)
while t<e    
ZW=0
h=maxx(P,f,n)
l=min(P,f,n)
Pcon=0.0
Pc=Center(P,h,n)     
Pr=reflection(P[h],Pc,P,n) 
if Check(Pr,P,l,f,n)   
            #f(Pr[1],Pr[2])<f(P[l][1],P[l][2]) 
    Pe=expansion(Pr,Pc,P,n)       
    if Check(Pe,P,h,f,n)   #f(Pe[1],Pe[2])<f(P[h][1],P[h][2]) 
            for i in 1:n-1
            P[h][i]=Pe[i]
            end           
    else        
            for i in 1:n-1
            P[h][i]=Pr[i]
            end           
    end          
elseif newmax(P,Pr,f,h,n)     
          if Check(Pr,P,h,f,n)   #f(Pr[1],Pr[2])<f(P[h][1],P[h][2])
          Ph=Pr          
          Pcon=contraction(Ph,Pc,n) 
          else       
            Pcon=contraction(P[h],Pc,n)
          end
            
          if Check(Pcon,P,h,f,n)   #f(Pcon[1],Pcon[2])<f(P[h][1],P[h][2])  
            for i in 1:n-1
            P[h][i]=Pcon[i]
            end      
          else 
            reduction(P,l,n)     
            ZW=1
          end   
else
     for i in 1:n-1
       P[h][i]=Pr[i]
      end         
end
#println(P)
if ZW==1
#P=Change(P)
#println(P)
end
e=count_e(P,n)
step+=1 
    if(step>100)
        e=t^2-1
    end
end

 for i in 1:n-1 
 print("Zmienna z")
 print(i)
 print(" = ")
 println(P[1][i])
    end
print("Steps = ")
print(step)
end







f(x,y,q,w,e,r,o,m,u,a,s,d)=x*x+y*y+q*q+w*w+e*e+o*o+u*u+a*a+s*s+d*d
m_Nelder(f,12)
