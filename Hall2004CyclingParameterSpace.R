#Autosomal suppression of meiotic drive can prevent sex chromosome cycling
#Code: Anjali Gupta
#Last worked on: 25 February 2024



#Set Working Directory
setwd("~/locate your directory")

# Parameters:
#   sdm = ssrm: Fitness cost of sex-ratio drive gene in males 
#   hd = hsr: Dominance of cost of drive in females
#   sdf = ssr: Fitness cost of sex-ratio drive gene in females
#   ha: Dominance of cost of autosomal suppressor
#   sa: Fitness cost of autosomal suppressor
#   sy: Fitness cost of Y-linked suppressor
#   d: Strength of drive
#   p1, p2, p3, ..., p9: Genotypic frequencies in females
#   q1, q2, q3, ..., q12: Genotypic frequencies in males
#   t: Start time (in generations)
#   t_final: End time (in generations)
#   Xdm = XSRm: Allelic frequency of driving X in males
#   Xdf = XSRf: Allelic frequency of driving X in females
#   Asm = Asupm: Allelic frequency of autosomal suppressor in males
#   Asf = Asupf: Allelic frequency of autosomal suppressor in females 
#   Ys = Ysup: Allelic frequency of Y-linked suppressor

#Checking if we get the same cycling behavior for the exact same parameters as Hall 2004 - region 4 (Figure 3)

#Making the dataframe for Hall 2004 - region 4 - Fig3 - Cycling parameter space

Hall_Region4 <- data.frame(ssrm=numeric(),
                           ssr=numeric(),
                           sy=numeric(),
                           hsr=numeric(),
                           d=numeric())

for(ssr in seq(0.1,0.4,0.1)) {
  for(sy in seq(0.1,0.5,0.1)) {
    
        ssrm = 0
        hsr = 3/5
        d = 2/5
    
        Hall_Region4 <- rbind(Hall_Region4, data.frame(ssrm=ssrm, ssr=ssr, sy=sy, hsr=hsr, d=d))
        
  }
}

for(ssr in seq(0.1,0.2,0.1)) {
  for(sy in seq(0.1,0.3,0.1)) {
        
        ssrm = 0
        hsr = 3/5
        d = 1/5
        
        Hall_Region4 <- rbind(Hall_Region4, data.frame(ssrm=ssrm, ssr=ssr, sy=sy, hsr=hsr, d=d))
        
  }
}

for(ssr in seq(0.1,0.3,0.1)) {
  for(sy in seq(0.1,0.7,0.1)) {
      
        ssrm = 0
        hsr = 2/5
        d = 2/5
        
        Hall_Region4 <- rbind(Hall_Region4, data.frame(ssrm=ssrm, ssr=ssr, sy=sy, hsr=hsr, d=d))
        
  }
}

for(ssr in seq(0.1,0.1,0.1)) {
  for(sy in seq(0.1,0.3,0.1)) {
        
        ssrm = 0
        hsr = 2/5
        d = 1/5
        
        Hall_Region4 <- rbind(Hall_Region4, data.frame(ssrm=ssrm, ssr=ssr, sy=sy, hsr=hsr, d=d))
        
  }
}

for(ssr in seq(0.1,0.2,0.1)) {
  for(sy in seq(0.4,0.7,0.1)) {
        
        ssrm = 0
        hsr = 0
        d = 2/5
        
        Hall_Region4 <- rbind(Hall_Region4, data.frame(ssrm=ssrm, ssr=ssr, sy=sy, hsr=hsr, d=d))
        
  }
}

# The while loop function simulates the dynamics of a population of a X-linked sex-ratio meiotic drive gene with Y-linked and autosomal suppressors over time.

# Parameters:
#   sdm = ssrm: Fitness cost of sex-ratio drive gene in males 
#   hd = hsr: Dominance of cost of drive in females
#   sdf = ssr: Fitness cost of sex-ratio drive gene in females
#   ha: Dominance of cost of autosomal suppressor
#   sa: Fitness cost of autosomal suppressor
#   sy: Fitness cost of Y-linked suppressor
#   d: Strength of drive
#   p1, p2, p3, ..., p9: Genotypic frequencies in females
#   q1, q2, q3, ..., q12: Genotypic frequencies in males
#   t: Start time (in generations)
#   t_final: End time (in generations)
#   Xdm = XSRm: Allelic frequency of driving X in males
#   Xdf = XSRf: Allelic frequency of driving X in females
#   Asm = Asupm: Allelic frequency of autosomal suppressor in males
#   Asf = Asupf: Allelic frequency of autosomal suppressor in females 
#   Ys = Ysup: Allelic frequency of Y-linked suppressor

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


library(ggplot2)

#Simulating our model over Hall's cycling parameter space
#Plots will be saved in the pdf

pdf("Hall2004CyclingParameterSpace.pdf")

for(i in 1:nrow(Hall_Region4)) {
  
  ssrm = Hall_Region4$ssrm[i]
  ssr = Hall_Region4$ssr[i]
  sy = Hall_Region4$sy[i]
  hsr = Hall_Region4$hsr[i]
  d = Hall_Region4$d[i]
  
  ha <- 0
  for(sa in 0) {
  
        p1=0.49 #XXAA ##a is sup ##A is wildtype 
        p2=0 #XXAa
        p3=0 #XXaa
        p4=0.01 #XsrXAA
        p5=0 #XsrXAa
        p6=0 #XsrXaa
        p7=0 #XsrXsrAA
        p8=0 #XsrXsrAa
        p9=0 #XsrXsraa
        q1=0.48 #XYAA
        q2=0 #XYAa
        q3=0 #XYaa
        q4=0.01 #XsrYAA
        q5=0 #XsrYAa
        q6=0 #XsrYaa
        q7=0 #XYsupAA
        q8=0 #XYsupAa
        q9=0 #XYsupaa
        q10=0.01 #XsrYsupAA
        q11=0 #XsrYsupAa
        q12=0 #XsrYsupaa
        
        t=0
        t_final=10000
        
        
        result_1 <- while_loop(ssrm,hsr,ssr,ha,sa,sy,d,
                               p1,p2,p3,p4,p5,p6,p7,p8,p9,
                               q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                               t,t_final)
        
        
        p1 <- ggplot(result_1, aes(x = t)) +
          geom_line(aes(y = Ysup), linetype = "solid", size = 0.5, color = "black") +
          geom_line(aes(y = XSRf), linetype = "dashed", size = 0.5, color = "blue") +
          geom_line(aes(y = XSRm), linetype = "solid", size = 0.5, color = "blue") +
          geom_line(aes(y = Asupm), linetype = "solid", size = 0.5, color = "green") +
          geom_line(aes(y = Asupf), linetype = "dashed", size = 0.5, color = "green") +
          xlim(0, 10000) +
          ylim(0, 1) +
          labs(x = "Time (generations)", y = "Frequency") +
          theme_bw() +
          ggtitle(paste("Asup (with cost) present ","ssrm",result_1$ssrm[1], "ssrf",result_1$ssr[1], "sa",result_1$sa[1], "sy",result_1$sy[1], "hsr",result_1$hsr[1], "ha",result_1$ha[1], "d",result_1$d[1])) +
          theme(axis.title = element_text(size=14),
                axis.text = element_text(size=8),
                title = element_text(size = 7)) +
          scale_color_manual(
            values = c("black", "blue", "blue", "green", "green"),
            labels = c("Ysup", "Xsr Female", "Xsr Male", "Asup Female", "Asup Male")
          ) +
          guides(
            color = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "solid", "dashed"))),
            linetype = guide_legend()
          )
        
        print(p1)
  }
}

dev.off()
