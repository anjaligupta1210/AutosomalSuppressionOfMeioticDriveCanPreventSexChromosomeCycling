#Simulation: Test invasion of Y-linked suppressor where the initial population is at equilibrium for X-linked driver and autosomal suppressor
#Code: Anjali Gupta
#Last worked on: 14 September 2023



# The while loop function simulates the dynamics of a population of a X-linked sex-ratio meiotic drive gene with Y-linked and autosomal suppressors over time.

# Parameters:
#   ssrm: Fitness cost of sex-ratio drive gene in males 
#   hsr: Dominance of cost of drive in females
#   ssr: Fitness cost of sex-ratio drive gene in females
#   ha: Dominance of cost of autosomal suppressor
#   sa: Fitness cost of autosomal suppressor
#   sy: Fitness cost of Y-linked suppressor
#   d: Strength of drive
#   p1, p2, p3, ..., p9: Genotypic frequencies in females
#   q1, q2, q3, ..., q12: Genotypic frequencies in males
#   t: Start time (in generations)
#   t_final: End time (in generations)
#   XSRm: Allelic frequency of driving X in males
#   XSRf: Allelic frequency of driving X in females
#   Asupm: Allelic frequency of autosomal suppressor in males
#   Asupf: Allelic frequency of autosomal suppressor in females 
#   Ysup: Allelic frequency of Y-linked suppressor

while_loop <- function(ssrm,hsr,ssr,ha,sa,sy,d,
                       p1,p2,p3,p4,p5,p6,p7,p8,p9,
                       q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                       t,t_final) {
  
  
  XSRm = (q10 + q11 + q12 + q4 + q5 + q6)/(q10 + q11 + q12 + q4 + q5 + q6 + q1 + q2 + q3 + q7 + q8 + q9) #freq of driver X in males
  XSRf = (p4 + p5 + p6 + 2*(p7 + p8 + p9))/(2*(p1 + p2 + p3) + p4 + p5 + p6 + p4 + p5 + p6 + 2*(p7 + p8 + p9)) #freq of driver X in females
  Asupm = (q11 + q2 + q5 + q8 + 2*(q12 + q3 + q6 + q9))/(q11 + q2 + q5 + 2*(q1 + q10 + q4 + q7) + q8 + q11 + q2 + q5 + q8 + 2*(q12 + q3 + q6 + q9)) #freq of Autosomal suppressor in males
  Asupf = (p2 + p5 + p8 + 2*(p3 + p6 + p9))/(p2 + p5 + 2*(p1 + p4 + p7) + p8 + p2 + p5 + p8 + 2*(p3 + p6 + p9)) #freq of Autosomal suppressor in females
  Ysup = (q10 + q11 + q12 + q7 + q8 + q9)/(q1 + q2 + q3 + q4 + q5 + q6 + q10 + q11 + q12 + q7 + q8 + q9) #freq of Y suppressor in males
  
  results = data.frame(ssrm=ssrm,hsr=hsr,ssr=ssr,ha=ha,sa=sa,sy=sy,d=d,
                       p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,
                       q1=q1,q2=q2,q3=q3,q4=q4,q5=q5,q6=q6,q7=q7,q8=q8,q9=q9,q10=q10,q11=q11,q12=q12,
                       t=t,
                       XSRm=XSRm,XSRf=XSRf,Asupm=Asupm,Asupf=Asupf,Ysup=Ysup)
  
  
  while(t<t_final &&
        !(q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12 == 0)) {
    
    u1 = v1 = 1
    u2 = v2 = 1 - ha*sa
    u3 = v3 = 1 - sa
    u4 = (1 - hsr*ssr)
    u5 = (1 - hsr*ssr)*(1 - ha*sa)
    u6 = (1 - hsr*ssr)*(1 - sa)
    u7 = (1 - ssr)
    u8 = (1 - ssr)*(1 - ha*sa)
    u9 = (1 - ssr)*(1 - sa)
    v4 = (1 - ssrm)
    v5 = (1 - ssrm)*(1 - ha*sa)
    v6 = (1 - ssrm)*(1 - sa)
    v7 = (1 - sy)
    v8 = (1 - sy)*(1 - ha*sa)
    v9 = (1 - sy)*(1 - sa)
    v10 = (1 - ssrm)*(1 - sy)
    v11 = (1 - ssrm)*(1 - sy)*(1 - ha*sa)
    v12 = (1 - ssrm)*(1 - sy)*(1 - sa)
    
    T = (v12*(p5/4 + p6/2 + p8/2 + p9)*(q8/4 + q9/2 + q11/4 + 
                                          q12/2)) + (v11*((p4/2 + p5/4 + p7 + p8/2)*(q8/4 + q9/2 + q11/4 +
                                                                                       q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q7/2 + q8/4 + q10/2 + 
                                                                                                                             q11/4))) + (v10*(p4/2 + p5/4 + p7 + p8/2)*(q7/2 + q8/4 + 
                                                                                                                                                                          q10/2 + q11/4)) + (v9*(p2/2 + p3 + p5/4 + p6/2)*(q8/4 + q9/2 + 
                                                                                                                                                                                                                             q11/4 + q12/
                                                                                                                                                                                                                             2)) + (v8*((p1 + p2/2 + p4/2 + p5/4)*(q8/4 + q9/2 + q11/4 + 
                                                                                                                                                                                                                                                                     q12/2) + (p2/2 + p3 + p5/4 + p6/2)*(q7/2 + q8/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                           q11/4))) + (v7*(p1 + p2/2 + p4/2 + p5/4)*(q7/2 + q8/4 + 
                                                                                                                                                                                                                                                                                                                                                       q10/2 + q11/4)) + (v6*(p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                          q5/4 + q6/
                                                                                                                                                                                                                                                                                                                                                                                                          2)) + (v5*((p5/4 + p6/2 + p8/2 + p9)*(q1/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  q2/4 + ((1/2) - d)*q4 + q5/4) + (p4/2 + p5/4 + p7 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     p8/2)*(q2/4 + q3/2 + q5/4 + q6/2))) + (v4*(p4/2 + p5/4 + p7 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  p8/2)*(q1/2 + q2/4 + ((1/2) - d)*q4 + q5/4)) + (v3*(p2/2 + p3 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        p5/4 + p6/2)*(q2/4 + q3/2 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q6/2)) + (v2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q6/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q2/4 + ((1/2) - d)*q4 + q5/4))) + (v1*(p1 + p2/2 + p4/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p5/4)*(q1/2 + q2/4 + ((1/2) - d)*q4 + q5/4)) + (u9*(p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       p6/2 + p8/2 + p9)*(q5/4 + q6/2 + q11/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            q12/2)) + (u8*((p4/2 + p5/4 + p7 + p8/2)*(q5/4 + q6/2 + q11/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(((1/2) + d)*q4 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q10/2 + q11/4))) + (u7*((p4/2 + p5/4 + p7 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         p8/2)*(((1/2) + d)*q4 + q5/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  q11/4))) + (u6*((p2/2 + p3 + p5/4 + p6/2)*(q5/4 + q6/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               q11/4 + q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             q8/4 + q9/2))) + (u5*((p1 + p2/2 + p4/2 + p5/4)*(q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                q6/2 + q11/4 + q12/2) + (p4/2 + p5/4 + p7 + p8/2)*(q2/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     q3/2 + q8/4 + q9/2) + (p2/2 + p3 + p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              p6/2)*(((1/2) + d)*q4 + q5/4 + q10/2 + q11/4) + (p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p6/2 + p8/2 + p9)*(q1/2 + q2/4 + q7/2 + q8/4))) + (u4*((p1 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           p2/2 + p4/2 + p5/4)*(((1/2) + d)*q4 + q5/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  q11/4) + (p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q8/4))) + (u3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q8/4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q9/2)) + (u2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q8/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q8/4))) + (u1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + q7/2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q8/4))
    
    q12next = (v12*(p5/4 + p6/2 + p8/2 + p9)*(q8/4 + q9/2 + q11/4 + 
                                                q12/2))/T
    q11next = (v11*((p4/2 + p5/4 + p7 + p8/2)*(q8/4 + q9/2 + q11/4 + 
                                                 q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q7/2 + q8/4 + q10/2 + 
                                                                                       q11/4)))/T
    q10next = (v10*(p4/2 + p5/4 + p7 + p8/2)*(q7/2 + q8/4 + q10/2 + 
                                                q11/4))/T
    q9next = (v9*(p2/2 + p3 + p5/4 + p6/2)*(q8/4 + q9/2 + q11/4 + q12/2))/T
    q8next = (v8*((p1 + p2/2 + p4/2 + p5/4)*(q8/4 + q9/2 + q11/4 + 
                                               q12/2) + (p2/2 + p3 + p5/4 + p6/2)*(q7/2 + q8/4 + q10/2 + 
                                                                                     q11/4)))/T
    q7next = (v7*(p1 + p2/2 + p4/2 + p5/4)*(q7/2 + q8/4 + q10/2 + q11/4))/T
    q6next = (v6*(p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + q5/4 + q6/2))/T
    q5next = (v5*((p5/4 + p6/2 + p8/2 + p9)*(q1/2 + 
                                               q2/4 + ((1/2) - d)*q4 + q5/4) + (p4/2 + p5/4 + p7 + 
                                                                                  p8/2)*(q2/4 + q3/2 + q5/4 + q6/2)))/T
    q4next = (v4*(p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + ((1/2) - d)*q4 +
                                              q5/4))/T
    q3next = (v3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q5/4 + q6/2))/T
    q2next = (v2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q5/4 + 
                                               q6/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + 
                                                                                    q2/4 + ((1/2) - d)*q4 + q5/4)))/T
    q1next = (v1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + ((1/2) - d)*q4 +
                                              q5/4))/T
    p9next = (u9*(p5/4 + p6/2 + p8/2 + p9)*(q5/4 + q6/2 + q11/4 + q12/2))/T
    p8next = (u8*((p4/2 + p5/4 + p7 + p8/2)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(((1/2) + d)*q4 + q5/4 + 
                                                                                     q10/2 + q11/4)))/T
    p7next = (u7*((p4/2 + p5/4 + p7 + p8/2)*(((1/2) + d)*q4 + q5/4 + 
                                               q10/2 + q11/4)))/T
    p6next = (u6*((p2/2 + p3 + p5/4 + p6/2)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + q8/4 + 
                                                                                     q9/2)))/T
    p5next = (u5*((p1 + p2/2 + p4/2 + p5/4)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p4/2 + p5/4 + p7 + p8/2)*(q2/4 + q3/2 + q8/4 + 
                                                                                     q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(((1/2) + d)*q4 + q5/4 + 
                                                                                                                          q10/2 + q11/4) + (p5/4 + p6/2 + p8/2 + p9)*(q1/2 + q2/4 + 
                                                                                                                                                                        q7/2 + q8/4)))/T
    p4next = (u4*((p1 + p2/2 + p4/2 + p5/4)*(((1/2) + d)*q4 + q5/4 + 
                                               q10/2 + q11/4) + (p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + 
                                                                                             q7/2 + q8/4)))/T
    p3next = (u3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q8/4 + q9/2))/T
    p2next = (u2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q8/4 + 
                                               q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                    q8/4)))/T
    p1next = (u1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + q7/2 + q8/4))/T
    
    
    
    
    
    XSRmnext = (q4next + q5next + q6next + q10next + q11next + q12next)/(q4next + q5next + q6next + q10next + q11next + q12next + q1next + q2next + q3next + q7next + q8next + q9next)
    XSRfnext = (p4next + p5next + p6next + 2*(p7next + p8next + p9next))/(p4next + p5next + p6next + 2*(p7next + p8next + p9next) + 2*(p1next + p2next + p3next) + p4next + p5next + p6next)
    Asupmnext = (q2next + q5next + q8next + q11next + 2*(q3next + q6next + q9next + q12next))/(q2next + q5next + q8next + q11next + 2*(q3next + q6next + q9next + q12next) + q2next + q5next + q8next + q11next + 2*(q1next + q4next + q7next + q10next))
    Asupfnext = (p2next + p5next + p8next + 2*(p3next + p6next + p9next))/(p2next + p5next + p8next + 2*(p3next + p6next + p9next) + p2next + p5next + p8next + 2*(p1next + p4next + p7next))
    Ysupnext = (q7next + q8next + q9next + q10next + q11next + q12next)/(q7next + q8next + q9next + q10next + q11next + q12next + q1next + q2next + q3next + q4next + q5next + q6next)
    
    t=t+1
    
    XSRm=XSRmnext
    XSRf=XSRfnext
    Ysup=Ysupnext
    Asupf=Asupfnext
    Asupm=Asupmnext
    
    p1=p1next
    p2=p2next
    p3=p3next
    p4=p4next
    p5=p5next
    p6=p6next
    p7=p7next
    p8=p8next
    p9=p9next
    q1=q1next
    q2=q2next
    q3=q3next
    q4=q4next
    q5=q5next
    q6=q6next
    q7=q7next
    q8=q8next
    q9=q9next
    q10=q10next
    q11=q11next
    q12=q12next
    
    results[length(results$t)+1,] <- c(ssrm=ssrm,
                                       hsr=hsr,
                                       ssr=ssr,
                                       ha=ha,
                                       sa=sa,
                                       sy=sy,
                                       d=d,
                                       p1=p1,
                                       p2=p2,
                                       p3=p3,
                                       p4=p4,
                                       p5=p5,
                                       p6=p6,
                                       p7=p7,
                                       p8=p8,
                                       p9=p9,
                                       q1=q1,
                                       q2=q2,
                                       q3=q3,
                                       q4=q4,
                                       q5=q5,
                                       q6=q6,
                                       q7=q7,
                                       q8=q8,
                                       q9=q9,
                                       q10=q10,
                                       q11=q11,
                                       q12=q12,
                                       t=t,
                                       XSRm=XSRm,
                                       XSRf=XSRf,
                                       Asupm=Asupm,
                                       Asupf=Asupf,
                                       Ysup=Ysup)
  }
  return(results)
}



result_all <- data.frame(ssrm=numeric(),hsr=numeric(),ssr=numeric(),ha=numeric(),sa=numeric(),sy=numeric(),d=numeric(),
                         p1mid=numeric(),p2mid=numeric(),p3mid=numeric(),p4mid=numeric(),p5mid=numeric(),p6mid=numeric(),p7mid=numeric(),p8mid=numeric(),p9mid=numeric(),
                         q1mid=numeric(),q2mid=numeric(),q3mid=numeric(),q4mid=numeric(),q5mid=numeric(),q6mid=numeric(),q7mid=numeric(),q8mid=numeric(),q9mid=numeric(),q10mid=numeric(),q11mid=numeric(),q12mid=numeric(),
                         XSRmmid=numeric(),XSRfmid=numeric(),Asupmmid=numeric(),Asupfmid=numeric(),Ysupmid=numeric(),
                         t=numeric(),
                         XSRm=numeric(),XSRf=numeric(),Asupm=numeric(),Asupf=numeric(),Ysup=numeric(),
                         p1=numeric(),p2=numeric(),p3=numeric(),p4=numeric(),p5=numeric(),p6=numeric(),p7=numeric(),p8=numeric(),p9=numeric(),
                         q1=numeric(),q2=numeric(),q3=numeric(),q4=numeric(),q5=numeric(),q6=numeric(),q7=numeric(),q8=numeric(),q9=numeric(),q10=numeric(),q11=numeric(),q12=numeric(),
                         DriveStable=character(),AutosomeStable=character(),YStable=character(),YInvade=character())


a <- seq(0,1,0.01)
h <- c(0,0.5,1)
b <- seq(0,0.5,0.01)

#This nested "for loop" loops over a range of possible parameters and tests for the invasion of a Y-linked suppressor in a population that is at equilibrium for X-linked drivers and autosomal suppressors

for(ssrm in 0){
  for(ssr in a){
    for(sa in c(0,0.1,0.2,0.5,0.9,1)){
      for(hsr in h){
        for(ha in h){
          for(sy in c(0.1,0.5,0.9)) {
            for(d in b) {
    
    p1=0.49 #XXAA ##a is autosomal suppressor ##A is standard autosome 
    p2=0 #XXAa
    p3=0 #XXaa
    p4=0.01 #XsrXAA
    p5=0 #XsrXAa ##Xsr is X-linked driver ##X is standard X 
    p6=0 #XsrXaa
    p7=0 #XsrXsrAA
    p8=0 #XsrXsrAa
    p9=0 #XsrXsraa
    q1=0.48 #XYAA
    q2=0 #XYAa
    q3=0 #XYaa
    q4=0.01 #XsrYAA
    q5=0.01 #XsrYAa
    q6=0 #XsrYaa
    q7=0 #XYsupAA ##Ysup is Y-linked suppressor ##Y is standard Y
    q8=0 #XYsupAa
    q9=0 #XYsupaa
    q10=0 #XsrYsupAA
    q11=0 #XsrYsupAa
    q12=0 #XsrYsupaa
    
    t=0
    t_final=5000
    
    result_1 <- while_loop(ssrm,hsr,ssr,ha,sa,sy,d,
                           p1,p2,p3,p4,p5,p6,p7,p8,p9,
                           q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                           t,t_final)
    
    
    p1=result_1$p1[5001]
    p2=result_1$p2[5001]
    p3=result_1$p3[5001]
    p4=result_1$p4[5001]
    p5=result_1$p5[5001]
    p6=result_1$p6[5001]
    p7=result_1$p7[5001]
    p8=result_1$p8[5001]
    p9=result_1$p9[5001]
    q1=result_1$q1[5001]
    q2=result_1$q2[5001]
    q3=result_1$q3[5001]
    q4=result_1$q4[5001]
    q5=result_1$q5[5001]
    q6=result_1$q6[5001]
    q7=result_1$q7[5001]
    q8=result_1$q8[5001]
    q9=result_1$q9[5001]
    q10=result_1$q10[5001]
    q11=result_1$q11[5001]
    q12=result_1$q12[5001]
    t=result_1$t[5001]
    XSRm=result_1$XSRm[5001]
    XSRf=result_1$XSRf[5001]
    Asupm=result_1$Asupm[5001]
    Asupf=result_1$Asupf[5001]
    Ysup=result_1$Ysup[5001]
    
    #This tests if the population has reached equilibrium and filters populations that have not gone extinct
    
    if((!is.na(t)==TRUE) &&
       t==5000){
      if(var(result_1$XSRm[(t-500):t])<1e-5 && 
         var(result_1$XSRf[(t-500):t])<1e-5 && 
         var(result_1$Asupm[(t-500):t])<1e-5 && 
         var(result_1$Asupf[(t-500):t])<1e-5 && 
         var(result_1$Ysup[(t-500):t])<1e-5) {
        
        p1mid=p1
        p2mid=p2
        p3mid=p3
        p4mid=p4
        p5mid=p5
        p6mid=p6
        p7mid=p7
        p8mid=p8
        p9mid=p9
        q1mid=q1
        q2mid=q2
        q3mid=q3
        q4mid=q4
        q5mid=q5
        q6mid=q6
        q7mid=q7
        q8mid=q8
        q9mid=q9
        q10mid=q10
        q11mid=q11
        q12mid=q12
        XSRmmid=XSRm
        XSRfmid=XSRf
        Asupmmid=Asupm
        Asupfmid=Asupf
        Ysupmid=Ysup
        
        #This loop tests for invasion of Y-linked suppressor and introduces the Y-linked suppressor into the population at equilibrium
        
        if (q4>=0.001) {
          q4 = q4 - 0.001
        }
        if (q4<0.001 & q3>=0.001) {
          q3 = q3 - 0.001
        }
        if (q4<0.001 & q3<0.001 & q2>=0.001) {
          q2 = q2 - 0.001
        }
        if (q4<0.001 & q3<0.001 & q2<0.001 & q1>=0.001) {
          q1 = q1 - 0.001
        }
        if (q4<0.001 & q3<0.001 & q2<0.001 & q1<0.001 & q5>=0.001) {
          q5 = q5 - 0.001
        }
        if (q4<0.001 & q3<0.001 & q2<0.001 & q1<0.001 & q5<0.001 & q6>=0.001) {
          q6 = q6 - 0.001
        }
        
        q10=0.001
        
        t=5001
        t_final=10000
        
        result_2 <- while_loop(ssrm,hsr,ssr,ha,sa,sy,d,
                               p1,p2,p3,p4,p5,p6,p7,p8,p9,
                               q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                               t,t_final)
        
        p1=result_2$p1[5000]
        p2=result_2$p2[5000]
        p3=result_2$p3[5000]
        p4=result_2$p4[5000]
        p5=result_2$p5[5000]
        p6=result_2$p6[5000]
        p7=result_2$p7[5000]
        p8=result_2$p8[5000]
        p9=result_2$p9[5000]
        q1=result_2$q1[5000]
        q2=result_2$q2[5000]
        q3=result_2$q3[5000]
        q4=result_2$q4[5000]
        q5=result_2$q5[5000]
        q6=result_2$q6[5000]
        q7=result_2$q7[5000]
        q8=result_2$q8[5000]
        q9=result_2$q9[5000]
        q10=result_2$q10[5000]
        q11=result_2$q11[5000]
        q12=result_2$q12[5000]
        t=result_2$t[5000]
        XSRm=result_2$XSRm[5000]
        XSRf=result_2$XSRf[5000]
        Asupm=result_2$Asupm[5000]
        Asupf=result_2$Asupf[5000]
        Ysup=result_2$Ysup[5000]
        
        
        if((!is.na(t)==TRUE) &&
           t==10000) { 
          
          #result_all dataframe contains information for populations and parameter spaces that we tested for invasion of the Y-linked suppressor
          
          result_all[length(result_all$ssrm)+1,]=c(ssrm,hsr,ssr,ha,sa,sy,d,
                                                   p1mid,p2mid,p3mid,p4mid,p5mid,p6mid,p7mid,p8mid,p9mid,
                                                   q1mid,q2mid,q3mid,q4mid,q5mid,q6mid,q7mid,q8mid,q9mid,q10mid,q11mid,q12mid,
                                                   XSRmmid,XSRfmid,Asupmmid,Asupfmid,Ysupmid,
                                                   t,XSRm,XSRf,Asupm,Asupf,Ysup,
                                                   p1,p2,p3,p4,p5,p6,p7,p8,p9,
                                                   q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                                                   DriveStable=if(var(result_2$XSRm[(5000-500):5000])<1e-5) "Yes" else "No",
                                                   AutosomeStable=if(var(result_2$Asupm[(5000-500):5000])<1e-5) "Yes" else "No",
                                                   YStable=if(var(result_2$Ysup[(5000-500):5000])<1e-5) "Yes" else "No",
                                                   YInvade=if(max(result_2$Ysup[3:length(result_2$Ysup)] > 0.005)) "Yes" else "No")
          
        
        }
      }
    }
            }
          }
        }
      }
    }
  }
}

#Write the output summmary file
write.csv(result_all,"YsupInvasion_AsupPresentInInitialPop_Fig1-2_S1-5.csv",row.names=FALSE)
