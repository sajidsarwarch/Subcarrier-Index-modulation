%Reference: (for alpha=1 only and SISO cases)   %%ref1: "Subcarrier Index Modulation for Spectral Efficient Frequency Division Multiplexing in Multi-Input Multi-Output Channels," 
%ref2: "Dual-Mode Index Modulation for Non-Orthogonal Frequency Division Multiplexing"
%Code is written by Muhammad Sajid Sarwar & I. Nyoman Apraz Ramatryana
clear;
clc;
%%%%%%%-----Define parameters----%%%%
G=32;% Number of groups
NF=128;%Total number of subcarriers
m=128;%Enter the total number of bits
CP=16;
N=NF/G;              %Number of subcarriers in each group
p=m/G;               %Number of bits per group
M=4;                 %M-ary Modulation
m_order=log2(M);
K=2;                 %Number of active subcarriers per group
L=10;                %channel tap


%%%%%-------- Signal generation ------%%%%%
v=1;
for SNRindB=0:2:20
    SNRindB
    if SNRindB<3
        Nn=40;               % data block
    elseif SNRindB<7
        Nn=400;
    elseif SNRindB<13
        Nn=1000;
    elseif SNRindB<17
        Nn=4000;
    else
        Nn=10000;
    end
    SNR=10^(SNRindB/10); %Signal-to-noise ratio
    sigma=sqrt(K*(NF+CP)/(2*NF*p*SNR)); % Noise variance
    B=Nn*m;
    num_err_bits(v)=0;
    
    for kk=1:Nn
        in_index=rand(1,(m*2+(m_order-1)*m/2))<0.5;
        signal=reshape(in_index,(2*p)+(m_order-1)*2,(G));%Every group of 4
        
            for g=1:G
                sig(:,g)=signal(:,g);%Divide the above-mentioned signals into two groups (corresponding to x_g: the unmodulated signal of the g-th block
                                     % (the number of bits in a group
                                     % should be according to
                                     % b=log2(C(N,K))+Klog2(Ma)+(N-K)log2(Mb)=b1+b2+b3
            end

        %%%%%%%%------Perform Mary modulation on b2 & b3--------%%%%%%
            for g=1:G
                p2_1(g)=bi2de(sig(3:3+m_order-1,g)');
                sym1(g)=pskmod(double(p2_1(g)),M); % The first Mary modulation symbol of the t group
                
            end  
        
        
            for g=1:G
                p2_2(g)=bi2de(sig(3+m_order:3+m_order+1,g)');             
                sym2(g)=pskmod(double(p2_2(g)),M);% The second Mary modulation symbol of the gth group
              
            end
      
        
        
            for g=1:G
                p2_11(g)=bi2de(sig(5+m_order:5+m_order+1,g)');
                
                sym11(g)=qammod(double(p2_11(g)),M);%  3rd modulation symbol of the gth group (could be a rotated version of psk. plz check ref2)
                
            end
      
        
        
       
            for g=1:G
                p2_22(g)=bi2de(sig(7+m_order:end,g)');
                
                sym22(g)=qammod(double(p2_22(g)),M);%  4th modulation symbol of the gth group (could be a rotated version of psk. plz check ref2)
                
            end

        
        %%%%%%%%----------IM modulation----------%%%%%
       
            for g=1:G
                
                C_temp= nchoosek ( 1 :N , K )  ; %generate number of possibe combinations C(N,K)
                i_C=bi2de([sig(2,g),sig(1,g)])+1;
                temp=zeros(N,1);
                temp(C_temp(i_C,1),1)=sym1(g);
                temp(C_temp(i_C,2),1)=sym2(g);
                
                %%%%% To allocate other Mode
                temp1=zeros(N,1);
                f=find(temp==0); % find out positions of unused subcarriers
                temp1(f(1,:),:)=sym11(g);
                temp1(f(2,:),:)=sym22(g);
%                 at=temp.*(1./temp);
%                 at(isnan(at))=0;
%                 temp1=zeros(N,1);
%                 for i=1:length(at)
%                     if at(i,:)==0
%                         temp1(i,:)=sym11(g,t);
%                         break;
%                     end
%                 end
%                 for j=length(at):-1:1
%                     if at(j,:)==0
%                         temp1(j,:)=sym22(g,t);
%                         break;
%                     end
%                 end
                
                %                 
                C0(:,g)=temp;% IM modulated symbols for  the G group (over nchoosek subcarriers)
                C1(:,g)=temp1; % IM modulated symbols for  the G group (over complementary nchoosek subcarriers)
            end

        C=C0+C1;
        %%%%%%%%-----Signal transmission and interleaving-----%%%%%
        for t=1:1
            x_F(:,:,t)=reshape((C(:,:,t).'),NF,1);%Transmit signal after interleaving
        end
        
        % Signal transmission
    
                
                h_T =  normrnd(0,sqrt(1/(2*2*L)),L,1)+1i*normrnd(0,sqrt(1/(2*2*L)),L,1); % ones(L,1);
                
                h_F =  fft([h_T;zeros(NF-L,1)]); %ones(128,1);

      
            Noise = 1*(normrnd(0,sigma,NF,1)+0*1i*normrnd(0,sigma,NF,1));
        
        
        %%%%%%%%%%%% signal reception %%%%%%%%%%%%%%
 
                yr_im_tmp= x_F .* h_F;
           
            yr_im=yr_im_tmp + Noise;%Received signal
        
        
        
        %Deinterleaving
        
            y_d= reshape((reshape(yr_im,G,N)).',NF,1);%Deinterleaved received signal y_d
       
        
 
                h_d = reshape((reshape(h_F,G,N)).',NF,1);% Channel coefficients after deinterleaving
 
        
        
            noise_d = reshape((reshape(Noise,G,N)).',NF,1);% Noise after deinterleaving
      
        
        %Grouping
        
        
            for g=1:G
                 y(:,:,1,g)=y_d((1+N*(g-1)):(N*g),:,1);
            end
      
        
       
                for g=1:G
                    h(:,:,1,g)=h_d(((1+N*(g-1)):(N*g)),:,1);% The following dimensions are equivalent to subscripts and subscripts
                end

        
        
            for g=1:G
                noise=noise_d((1+N*(g-1)):(N*g),:);
            end
      
     
        %% ------------%% Joint Signal Detection
        
        % Create a lookup table (transmission symbols obtained after IM mapping and Mary modulation of the original bit information)
        lut=OFDM_DM_Look_Up_Table(N,K,M);
        
        for g=1:G
            for a=1:length(lut) % lookup table column
             
                    tmp1=diag(lut(:,a))*h(:,:,1,g);  %% The received signal joint detection
                    met1(a)=(norm(y(:,:,1,g)-tmp1,'fro')).^2;

                    
                    met(a)=met1(a); 

                end
          
            [out_T1]=find(met == min(min(met)));%apply maximum likelihood detection (joint decoding of index and symbols)
            
           %%%%%%%-----------bits conversion ---------------%%%%%%%%%%
            temp_A=(de2bi((out_T1-1),8+(m_order-1)*2,'left-msb')).'; %%%%% the no. of bits (b1+b2+b3) should be equal to the bits transmitted by a group. If not, adjust accordingly.
            temp_A([3:3+m_order-1])=temp_A(flipud([3:3+m_order-1]'));
            temp_A([3+m_order:3+m_order+m_order-1])=temp_A(flipud([3+m_order:3+m_order+m_order-1]'));
            temp_A([3+m_order+m_order:3+m_order+m_order+m_order-1])=temp_A(flipud([3+m_order+m_order:3+m_order+m_order+m_order-1]'));
            temp_A([3+m_order+m_order+m_order:3+m_order+m_order+m_order+m_order-1])=temp_A(flipud([3+m_order+m_order+m_order:3+m_order+m_order+m_order+m_order-1]'));
           
            out(:,g) = temp_A;% Corresponds to the g-th group of bits 
           
            
            
        end
        out_index=reshape(out,1,[]);
     
        out_index=[out_index];
        
        
        %% ------Calculate the number of bit errors-----%%
        for num=1:(m)
            if in_index(num) ~= out_index(num)
                num_err_bits(v) = num_err_bits(v) + 1;
            end
        end
        kk
    end
    prb(v) = num_err_bits(v)/B;
    v=v+1;
end
SNRindB=0:2:20;
semilogy(SNRindB,prb,'->k');