clc;
clearvars;
% ---------------------------ENCODING------------------------------%
inp_len=input('enter number of input bits');
%inp_len=500;
K=input('enter constraint length');
%K=3;
zer=zeros(1,K-1);                      % z is the matrix of zeros to be added

for i=1:inp_len
in(i)=rand <0.5;
end

enc_in=[zer,in,zer]; % enc_in is the encoded input with zeros
n=input('enter number of parity bits per window');
%n=2;
gen=zeros(n,K); % gen is a matrix which stores the generator polynomials

for i=1:n
g=input('enter generator polynomial','s');
o=char(num2cell(g));
o=reshape(str2num(o),1,[]);
       
         for j=1:K
         gen(i,j)=o(1,j);
         end
end

enc_bits=zeros(1,n*(inp_len+K-1));

for i=1:inp_len+K-1
    for j=1:n
        for d=K:-1:1
            enc_bits(1,j+(i-1)*n)=enc_in(1,K-d+i)*gen(j,d) + enc_bits(1,j+(i-1)*n);
        end
        enc_bits(1,n*(i-1)+j)=mod(enc_bits(1,n*(i-1)+j),2);
    end
end

display(enc_bits); % enc_bits are the bits which are encoded
%-------------------------CHANNEL ENCODING----------------------------
gamaset=0:0.5:8;
bsc_err=zeros(size(gamaset));
bec_err=zeros(size(gamaset));
gau_err=zeros(size(gamaset));
final=zeros(1,n*(K+inp_len-1));

for chann=1:3
      vr=0;
         for gama=gamaset
             vr=vr+1;
             gamalinear=10.^(gama./10);
             pr(vr)=qfunc(sqrt(gamalinear*2/n));
             sigma=sqrt(n/2.*gamalinear);
%--------------------Stimulation---------------------%             
             for ksim=1:2000
             if(chann==1) % for BSC
                   for i=1:n*(inp_len+K-1)
                          if(rand<pr(vr))
                             if(enc_bits(1,i)==1)
                                final(1,i)=0;
                             else
                                final(1,i)=1;
                             end
                          else
                                final(1,i)=enc_bits(1,i);
                          end
                   end
             end
             
             if(chann==2) % for BEC
              for i=1:n*(K+inp_len-1)
                  if(rand<pr(vr))
                  final(1,i)=-300;
                  else
                  final(1,i)=enc_bits(1,i);
                  end
              end
             end
             
            if(chann==3) % for Gaussian
             sigma=sqrt(n/(2*gamalinear));
              
                  final=enc_bits;
                  noise=sigma*randn(1,length(final));
                   sign=2*final-1;
            
                 rec=sign + noise;
                 final=rec;
           end
   

    sou=final;
%--------------------decoder------------%    
if(chann==3)
    h=-1;
else
    h=0;
end

Box=zeros(4,n*(inp_len+K-1)/2+1); % Box is a matrix which stores the hamming values of the nodes
diff_lin=zeros(8,n*(inp_len+K-1)/2); % diff_lin is a matrix which stores the hamming distance of all branches
bits=[h,h;1,1;1,h;h,1;1,1;h,h;h,1;1,h]; % all bits as per their allocated branch number
 
      
         for i=1:2:n*(inp_len+K-1)-1
             for j=1:8
              diff_lin(j,(i+1)/2)=abs(sou(i)-bits(j,1));
              diff_lin(j,(i+1)/2)=abs(sou(i+1)-bits(j,2))+diff_lin(j,(i+1)/2);
             end
         end
     

if(chann==3)
    diff_lin=diff_lin.*diff_lin;
end
Box(2,1)=[100000];        % randomly assigned the starting unused boxes
Box(3,1)=[100000];
Box(4,1)=[100000];
Box(2,3)=[100000];
Box(2,4)=[100000];

rd_bits=zeros(1,n*(inp_len+K-1)/2); % matrix which stores the values of the decoded bits in reverse order

for i=2:n*(inp_len+K-1)/2+1
      for j=1:4
          if(j==1)
              if((Box(1,i-1)+diff_lin(1,i-1))<=(Box(3,i-1)+diff_lin(5,i-1)))
                    Box(1,i)=(Box(1,i-1)+diff_lin(1,i-1));
              else
                    Box(1,i)=(Box(3,i-1)+diff_lin(5,i-1));
              end
          end

          if(j==2)
              if((Box(1,i-1)+diff_lin(2,i-1))<=(Box(3,i-1)+diff_lin(6,i-1)))
                    Box(2,i)=(Box(1,i-1)+diff_lin(2,i-1));
              else
                    Box(2,i)=(Box(3,i-1)+diff_lin(6,i-1));
              end
          end
      
          if(j==3)
             if((Box(2,i-1)+diff_lin(3,i-1))<=(Box(4,i-1)+diff_lin(7,i-1)))
                  Box(3,i)=(Box(2,i-1)+diff_lin(3,i-1));
             else
                  Box(3,i)=(Box(4,i-1)+diff_lin(7,i-1));
             end
          end
         if(j==4)
            if((Box(2,i-1)+diff_lin(4,i-1))<=(Box(4,i-1)+diff_lin(8,i-1)))
                  Box(4,i)=(Box(2,i-1)+diff_lin(4,i-1));
            else
                  Box(4,i)=(Box(4,i-1)+diff_lin(8,i-1));
            end
         end
     end
end

  if(Box(1,n*(inp_len+K-1)/2+1)<=Box(2,n*(inp_len+K-1)/2+1) && Box(1,n*(inp_len+K-1)/2+1)<=Box(3,n*(inp_len+K-1)/2+1) && Box(1,n*(inp_len+K-1)/2+1)<=Box(4,n*(inp_len+K-1)/2+1))
  j=1;
  elseif(Box(2,n*(inp_len+K-1)/2+1)<Box(1,n*(inp_len+K-1)/2+1) && Box(2,n*(inp_len+K-1)/2+1)<=Box(3,n*(inp_len+K-1)/2+1) && Box(2,n*(inp_len+K-1)/2+1)<=Box(4,n*(inp_len+K-1)/2+1))
  j=2;
  elseif(Box(3,n*(inp_len+K-1)/2+1)<Box(2,n*(inp_len+K-1)/2+1) && Box(3,n*(inp_len+K-1)/2+1)<Box(1,n*(inp_len+K-1)/2+1) && Box(3,n*(inp_len+K-1)/2+1)<=Box(4,n*(inp_len+K-1)/2+1))
  j=3;
  else
  j=4;
  end
  p=0;

  for i=n*(inp_len+K-1)/2+1:-1:2
  p=p+1;
  if(j==1)
      if(Box(1,i)==Box(1,i-1)+diff_lin(1,i-1))
          j=1;
          rd_bits(1,p)=0;
      else
         j=3;
         rd_bits(1,p)=0;
      end
 elseif(j==2)
     if(Box(2,i)==Box(1,i-1)+diff_lin(2,i-1))
        j=1;
        rd_bits(1,p)=1;
     else
       j=3;
       rd_bits(1,p)=1;
     end
 elseif(j==3)
    if(Box(3,i)==Box(2,i-1)+diff_lin(3,i-1))
       j=2;
       rd_bits(1,p)=0;
    else
      j=4;
      rd_bits(1,p)=0;
    end   
    
  elseif(j==4)
        if(Box(4,i)==Box(2,i-1)+diff_lin(4,i-1))
            j=2;
            rd_bits(1,p)=1;
        else
            j=4;
            rd_bits(1,p)=1;
        end
      end
end

for i=1:inp_len
    dec_bits(i)=rd_bits(inp_len+K-i);
end

for nb=1:inp_len
    if(chann==1)
        if(in(nb)~=dec_bits(nb))
        bsc_err(vr)=bsc_err(vr)+1;
        end
    end
    if(chann==2)
        if(in(nb)~=dec_bits(nb))
        bec_err(vr)=bec_err(vr)+2.4;
        end
    end
    if(chann==3)
        if(in(nb)~=dec_bits(nb))
        gau_err(vr)=gau_err(vr)+1;
        end
    end
    end
             end
    end

    p_err_gau=gau_err/(2000*inp_len);
    p_err_bsc=bsc_err/(2000*inp_len);
    p_err_bec=bec_err/(2000*inp_len);
    
    semilogy(gamaset,p_err_bsc,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
    hold on;
    semilogy(gamaset,p_err_gau,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
    hold on;
    semilogy(gamaset,p_err_bec,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
    xlabel('SNR per Bit in dB'); ylabel('Probability of Bit Error'); grid
    legend('BSC','Gaussian Noise','BEC');     

end
