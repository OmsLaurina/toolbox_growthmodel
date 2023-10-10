def growth_model_2P1Z_v10_simplify(Psupply,time,**kwargs):
    
    global dt, umax1, umax2, gmax1, gmax2, kP1, kP2, kZ, eZ, epsilon, mZ, gamma1, gamma2, u1, u2, g1, g2, PP1, PP2, G1, G2, P1, P2, Z, PO4, P01, P02, Z0, PO40

    """" definition of model variables """
    
    ##Creating empty vectors
    
    # growth rates [d^-1]
    u1=[]
    u2=[]		
    # grazing rates [d^-1]
    g1=[]
    g2=[]	
			
    #Primary production [mmolC^{-1} m^{3} d^{-1}]
    PP1=[]
    PP2=[]
    
    #grazing [mmolC^{-1} m^{3} d^{-1}]
    G1=[]
    G2=[]
    
    #excretion [mmolC^{-1} m^{3} d^{-1}]
    exc=[]
    
    #death [mmolC^{-1} m^{3} d^{-1}]
    #Natural mortality
    d_P1 = []
    d_P2 = []
    d1_Z=[]
    #High trophic level mortality 
    d2_Z = []
    
    #recycling
    rec_P1 = []
    rec_P2 = []
    rec_Z = []
    rec_exc = []
    
    #Washout (biomass out of the system) [mmolC^{-1} m^{3} d^{-1}]
    w_P1 = []
    w_P2 = []
    w1_Z = []
    w2_Z = []
    
    Export=[]
    
    # biomass of Phytoplankton, Zooplankton and PO4 [mmolC/m³]
    P1=[] 			
    P2=[] 						
    Z = []
    PO4 = []
    
    """ default value of model parameters and initial values """
    
  
    defaultKwargs = {
        'umax1':1.9872,   # maximum growth rates of P1 [d^{-1}]
        'umax2':2.7648,   # maximum growth rates of P2 [d^{-1}]
        'gmax1':3.89,     # maximum grazing rates of Z on P1 [d^{-1}]
        'gmax2':0.43,     # maximum grazing rates of Z on P2 [d^{-1}]
        'kP1':1,          # half-saturation constant for P1 on PO4 [mmolC m^{-3}]
        'kP2':3,	      # half-saturation constant for P2 on PO4 [mmolC m^{-3}]
        'kZ1':5,		  # half-saturation constant for Z on P1 [mmolC m^{-3}]
        'kZ2':20,		  # half-saturation constant for Z on P2 [mmolC m^{-3}]
        'mP1':0.1,	      # P1 natural mortality rate [d^{-1}]
        'mP2':0.2,	      # P2 natural mortality rate [d^{-1}]
        'm1Z':0.1,	      # Z natural mortality rate [d^{-1}]
        'm2Z':0.061,	  # Z quadratic mortality rate [mmolC^{-1} m^{3} d^{-1}]
        'gamma':0.6,      # conversion factor from P to Z [/]
        'epsilonP':1,     # fraction of P natural mortality that is available as regenerated PO4 [/]
        'epsilon1Z':0.3,  # fraction of Z natural mortality that is available as regenerated PO4 [/]
        'epsilon2Z':0.7,  # fraction of Z excretion that is available as regenerated PO4 [/]
        'P1_ini':0.6,     # initial biomass of P1 [mmolC m^{-3}]
        'P2_ini':0.1,     # initial biomass of P2 [mmolC m^{-3}]
        'Z_ini':0.6,      # initial biomass of Z [mmolC m^{-3}]
        'PO4_ini':0.5,    # initial biomass of PO4 [mmolC m^{-3}]
        'DT':0.2          # timestep of discretisation  [day]
        }
    
    arg = { **defaultKwargs, **kwargs}
    
    dt=arg['DT']
    
    # Initial conditions at time=0			
    P1.append(arg['P1_ini'])
    P2.append(arg['P2_ini'])
    Z.append(arg['Z_ini'])
    PO4.append(arg['PO4_ini'])
    
    
    nb_time=len(time)
    
    for t in range(0,nb_time-1):

        # growth rates (monod function)
        u1next=PO4[t-1]/(arg['kP1'])*arg['umax1']
        u2next=PO4[t-1]/(arg['kP2'])*arg['umax2']
       
        u1.append(u1next)
        u2.append(u2next)
        
        # functional response = grazing rates (holling type II)
        g1next=P1[t-1]/(arg['kZ1'])*arg['gmax1']
        g2next=P2[t-1]/(arg['kZ2'])*arg['gmax2']
       
        g1.append(g1next)
        g2.append(g2next)
   
   	    # FLUXES
          
        # UPTAKE (PRIMARY PRODUCTION)
        PP1next=u1[t]*P1[t-1]
        PP2next=u2[t]*P2[t-1]
       
        PP1.append(PP1next)
        PP2.append(PP2next)
       
        # GRAZING
        G1_next=g1[t]*Z[t-1]
        G2_next=g2[t]*Z[t-1]
        
        G1.append(G1_next)
        G2.append(G2_next)
        
        # EXCRETION
        excnext=(1-arg['gamma'])*g1[t]*Z[t-1]+(1-arg['gamma'])*g2[t]*Z[t-1]
        
        exc.append(excnext)
        
        # MORTALITY
        d_P1next = arg['mP1']*P1[t-1]
        d_P2next = arg['mP2']*P2[t-1]
        d1_Znext = arg['m1Z']*Z[t-1]
        d2_Znext = arg['m2Z']*Z[t-1]**2
        
        d_P1.append(d_P1next)
        d_P2.append(d_P2next)
        d1_Z.append(d1_Znext)
        d2_Z.append(d2_Znext)
        
        # RECYCLING
        rec_P1next = arg['epsilonP']*d_P1[t]
        rec_P2next = arg['epsilonP']*d_P2[t]
        rec_Znext = arg['epsilon1Z']*d1_Z[t]
        rec_excnext = arg['epsilon2Z']*exc[t]
        
        rec_P1.append(rec_P1next)
        rec_P2.append(rec_P2next)
        rec_Z.append(rec_Znext)
        rec_exc.append(rec_excnext)
        
        #WASHOUT
        #eps = Psupply[t]/(arg['m1Z']*Z[t-1])
        
        w_P1next = (1-arg['epsilonP'])*d_P1[t]
        w_P2next = (1-arg['epsilonP'])*d_P2[t]
        w1_Znext = (1-arg['epsilon1Z'])*d1_Z[t]
        w2_Znext = (1-arg['epsilon2Z'])*exc[t]
        
        w_P1.append(w_P1next)
        w_P2.append(w_P2next)
        w1_Z.append(w1_Znext)
        w2_Z.append(w2_Znext)
        
        Exportnext = w_P1[t]+w_P2[t]+w1_Z[t]+w2_Z[t]+d2_Z[t]
        
        Export.append(Exportnext)
    
        #P1 puise en premier dans le pool de nutriments
        max_available_PO4 = PO4[t-1]+Psupply[t]*dt+rec_P1[t]*dt+rec_P2[t]*dt+rec_Z[t]*dt+rec_exc[t]*dt
       
        if PP1[t]>max_available_PO4/dt:
            PP1[t]=max_available_PO4/dt
            
        max_available_PO4=max_available_PO4-PP1[t]*dt
        if PP2[t]>max_available_PO4/dt:
            PP2[t]=max_available_PO4/dt  
        
        # BIOMASS
        P1_next=P1[t-1]+PP1[t]*dt-G1[t]*dt-d_P1[t]*dt
        if P1_next<=0:
            P1_next=0
        
        P2_next=P2[t-1]+PP2[t]*dt-G2[t]*dt-d_P2[t]*dt
        if P2_next<=0:
            P2_next=0
            
        Z_next=Z[t-1]+arg['gamma']*G1[t]*dt+arg['gamma']*G2[t]*dt-d1_Z[t]*dt-d2_Z[t]*dt
        if Z_next<=0:
            Z_next=0
            
        PO4_next=PO4[t-1]+Psupply[t]*dt+rec_P1[t]*dt+rec_P2[t]*dt+rec_Z[t]*dt+rec_exc[t]*dt-PP1[t]*dt-PP2[t]*dt
        if PO4_next<=0:
            PO4_next=0
        
        P1.append(P1_next)
        P2.append(P2_next)
        Z.append(Z_next)
        PO4.append(PO4_next)
        
        
    return P1,P2,Z,PO4,Export,u1,u2,g1,g2, arg