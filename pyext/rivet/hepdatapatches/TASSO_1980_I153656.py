import yoda,math
# fix bin widths and transform axis to x for relevant dists (x02)
def patch(path, ao):
    if "TASSO_1980_I153656" in path :
        temp=path.split("/")[-1].split("-")
        i = int(temp[0][1:])
        if i<2 or i>13 : return ao
        j = int(temp[2][1:])
        if (i>1 and i<8) and j>2 : return ao
        if(i>=2 and i<=4) : sqrtS=12.
        else              : sqrtS=30.
        for k in range(0,len(ao.points())) :
            step2=-1
            if(i==2) :
                if(k<3) :
                    step1 = 0.05
                elif(k<4) :
                    step1 = 0.2
                else :
                    step1 = 0.3
            elif(i==3) :
                if(k<2) :
                    step1 = 0.05
                else :
                    step1 = 0.2
            elif(i==4) :
                if(k<2) :
                    step1 = 0.125
                else :
                    step1 = 0.6
            elif(i==5) :
                if(k>4) :
                    step1 = 0.1
                elif(k==4) :
                    step2=0.1
                    step1=0.025
                elif(k==3) :
                    step2=0.025
                    step1=0.05
                elif(k==2) :
                    step1=0.025
                    step2=0.05
                else :
                    step1=0.025
            elif(i==6) :
                if(k>2) :
                    step1 = 0.1
                elif(k==2) :
                    step1=0.025
                    step2=0.1
                elif(k==1) :
                    step2=0.025
                    step1=0.05
                else :
                    step1=0.05
            elif(i==7) :
                if(k<2) :
                    step1=0.1
                else :
                    step1=1.665-0.9
            elif(i==8 or i==9) :
                if(k<2) :
                    step1=0.05
                else :
                    step1=0.2
            elif(i==10 or i==13) :
                if(k==0) :
                    step1=0.125
                elif(k==1) :
                    step1=0.125
                    step2=.15
                else :
                    step1=ao.points()[k].xErrs()[0]
                    step2=ao.points()[k].xErrs()[1]
            elif(i==11 or i==12) :
                if(k==0) :
                    step1=0.05
                elif(k==1) :
                    step1=0.05
                    step2=0.025
                elif(k==2) :
                    step1=0.025
                    step2=0.1
                else :
                    step1=0.1
            if(step2<0) : step2 = step1
            ao.points()[k].setXErrs((step1,step2))
            if(j==2) :
                 ao.points()[k].setXErrs((2.*step1/sqrtS,2.*step2/sqrtS))
                 ao.points()[k].setX( ao.points()[k].x()*2./sqrtS)


    return ao
