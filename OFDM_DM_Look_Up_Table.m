function [lut]=OFDM_DM_Look_Up_Table(N,K,M)
% clc, clear all
% N=4; M=4;K=2;
C= nchoosek ( 1 :N , K )  ;
symbol=pskmod(0:M-1,M);
symbol_1=qammod(0:M-1,M);
lut=[];
for i_C=1:N
    for i_M1=1:length(symbol)
        for i_M2=1:length(symbol)
            for i_M3=1:length(symbol)
                for i_M4=1:length(symbol)
                    temp=zeros(N,1);
                    temp(C(i_C,1),1)=symbol(i_M1);
                    temp(C(i_C,2),1)=symbol(i_M2);
                    
                    f=find(temp==0);   %%%%% to find the unused subcarriers after 1st mode
                    %         for i=1:length(f)
                    temp(f(1,:),:)=symbol_1(i_M3);
                    temp(f(2,:),:)=symbol_1(i_M4);
                    lut=[lut,temp];
                    
                    
                end
            end
        end
    end
end
end