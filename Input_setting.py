

def basic_sys_setting(Choose):

    if Choose == 1:
        Nk = 3
        #m = 50
        #m = 5
        #m = 10
        #m = 20
        #m = 40
        m = 64

        #mea = 100                   #  mea =  number of shots
        #mea = 500                   #  mea =  number of shots    
        mea = 2400                   #  mea =  number of shots    
            
    elif Choose == 2:
        Nk = 4
        #m = 50
        #m = 100
        #m = 200
        m = 256

        #mea = 400                   #  mea =  number of shots
        #mea = 800                   #  mea =  number of shots
        #mea = 1200                   #  mea =  number of shots
        #mea = 2400                   #  mea =  number of shots

        mea = 4800                   #  mea =  number of shots

        #mea = 8600                   #  mea =  number of shots

    elif Choose == 3:
        Nk = 6
        #m = 50
        #m = int(0.05* (4**Nk))        # 204
        #m = int(0.1* (4**Nk))        #  409
        m = int(0.2* (4**Nk))       #  819
        #m = int(0.3* (4**Nk))       #  1228         
        #m = int(0.4* (4**Nk))       #  1638
        #m = int(0.5* (4**Nk))       #  2048
        #m = int(0.6* (4**Nk))       #  2457

        #mea = 1200                   #  mea =  number of shots
        #mea = 2400                   #  mea =  number of shots
        #mea = 3600                   #  mea =  number of shots
        #mea = 4800                   #  mea =  number of shots
        mea = 8600                    #  mea =  number of shots

        #mea = 9600                   #  mea =  number of shots
        #mea = 19200                   #  mea =  number of shots

    elif Choose == 35:
        Nk = 7
        #m = 50
        #m = int(0.05* (4**Nk))        #    819
        #m = int(0.1* (4**Nk))        #    1638
        m = int(0.2* (4**Nk))          #   3276 

        mea = 8600                    #  mea =  number of shots

        #mea = 9600                   #  mea =  number of shots

    elif Choose == 4:
        Nk = 8
        #m = 50
        #m = int(0.05* (4**Nk))      # 3276
        m = int(0.1* (4**Nk))        # 6553
        #m = int(0.2* (4**Nk))         # 13107
        #m = int(0.3* (4**Nk))         # 19660
        #m = int(0.4* (4**Nk))

        #mea = 1200                   #  mea =  number of shots
        #mea = 2400
        #mea = 4800
        #mea = 6000
        #mea = 6400
        #mea = 8000
        mea = 8600

    elif Choose == 5:
        Nk = 10
        #m = 50
        #m = int(0.05* (4**Nk))        #   52428
        #m = int(0.1* (4**Nk))        #  104857
        #m = int(0.2* (4**Nk))        #  209715
        #m = int(0.25* (4**Nk))        #  262144
        m = int(0.3* (4**Nk))        #  314572


        #mea = 1200                   #  mea =  number of shots
        mea = 8600

    elif Choose == 7:
        Nk = 12
        m = int(0.05* (4**Nk))        # 838860
        #m = int(0.1* (4**Nk))          # 1677721    

        mea = 8600
        #mea = 17200

    elif Choose == 8:
        Nk = 13
        m = int(0.05* (4**Nk))       
        #m = int(0.1* (4**Nk))       

        mea = 8600

    return Nk, m, mea
