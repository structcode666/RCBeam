from dataclasses import dataclass
from typing import Optional, List
import math



def load_combos(pl_dead : float, pl_live : float, udl_dead : float, udl_live:float, conc_beam = None) -> []:

    ##Load_combos##
    pl_short_term = pl_dead + 0.7*pl_live
    pl_long_term = pl_dead + 0.4*pl_live
    pl_ultimate = 1.2*pl_dead + 1.5*pl_live

    udl_short_term = udl_dead + conc_beam.get_self_weight() + 0.7*udl_live
    udl_long_term = udl_dead + conc_beam.get_self_weight() + 0.4*udl_live
    udl_ultimate = 1.2*(udl_dead+conc_beam.get_self_weight()) + 1.5*udl_live

    load_combo_list = [pl_ultimate, 
                       udl_ultimate, 
                       pl_short_term, 
                       pl_long_term, 
                       udl_short_term, 
                       udl_long_term]

    return load_combo_list



def beam_analysis(span : float, load_combo_list = []):

    ## Get Bending Moment##
    ult_des_bending_moment = (load_combo_list[0]*span)/4 + (load_combo_list[1]*span*span)/8
    ult_shear = load_combo_list[0]/2 + (load_combo_list[1]*span*0.5)
    service_des_bending_moment = (load_combo_list[2]*span)/4 + (load_combo_list[4]*span*span)/8
    service_des_bending_moment_long_term = (load_combo_list[3]*span)/4 + (load_combo_list[5]*span*span)/8

    return ult_des_bending_moment, ult_shear, service_des_bending_moment


@dataclass
class concrete_rect_beam:
    """
    Dataclass that represents a rectangular concrete beam with top and bottom reinforcement
    """
    b : float
    d : float
    conc_grade : int 
    des_bending_moment : Optional[float] = None
    ser_design_moment : Optional[float] = None
    ser_design_bending_moment_long_term : Optional[float]= None
    f_sy: Optional[int] = 500.
    As_top : Optional[float] = 0.
    steel_modulus = 200000.
    
    
    def get_self_weight(self):

        self_weight = (self.b*self.d*25)/(1000*1000)

        return self_weight
   
    def get_Ig(self):
        """
        Gives concrete rectangular section gross properties
        """
        conc_sec_Ig = self.b*self.d*self.d*self.d/12
        return conc_sec_Ig
    
    def get_Ief(self):

         if (self.ult_tensile_reo() == 'Increase beam depth'):
            return "Increase Beam Depth"
         
         else:
            #web reo ratio for tensile reo#
            reo_ratio_tension = (self.ult_tensile_reo())/(self.b*self.d)  #Note web reo ratio for compression reo has been ignored for SS beam - can improve for continous beams##

            ## Ief##
            if reo_ratio_tension >= 0.001*(self.conc_grade)**(1/3):
                Ief = min((self.b*(self.d)**3*((5-0.04*self.conc_grade)*reo_ratio_tension+0.002)) , 0.1*(self.b)*(self.d)**3)
            elif reo_ratio_tension < 0.001*(self.conc_grade)**(1/3):
                Ief = min(self.b*(self.d)**3*(0.055*(self.conc_grade)**(1/3) - 50*reo_ratio_tension),0.06*self.b*(self.d)**3 )

            return Ief
    
    def get_kcs(self):

        kcs = max(2-1.2*(self.As_top/concrete_rect_beam.ult_tensile_reo(self)), 0.8)

        return kcs
    
    def get_Ec(self):
        
        if self.conc_grade == 20:
            Ec = 24000
        elif self.conc_grade == 25:
            Ec = 26700
        elif self.conc_grade == 32:
            Ec = 30100
        elif self.conc_grade == 40:
            Ec = 32800
        elif self.conc_grade == 50:
            Ec = 34800
        elif self.conc_grade == 65:
            Ec = 37500
        elif self.conc_grade == 80:
            Ec = 39600
        
        return Ec
        
    
    def ult_tensile_reo(self):
        
        #Caluclate stress block modification factors#
        alpha_2 = max(0.67, 0.85-0.0015*self.conc_grade)
        gamma = max(0.67, 0.97-0.0025*self.conc_grade)
        
        #Check assumed Astd#
        lever_arm_assumed = 0.9 * self.d
        Ast_assumed = (self.des_bending_moment*1000*1000)/(0.85*500*lever_arm_assumed)
    
        #Calculate ku#
        k_u = (Ast_assumed*self.f_sy)/(alpha_2 * self.conc_grade * gamma * self.d * self.b)

        if k_u>0.36:
            return "Increase beam depth"
        elif k_u<=0.36:
            #Calculate final lever arm
            lever_arm_final = self.d - (0.5*gamma*k_u*self.d)
            
            #Capacity Reduction Factor#
            phi = min((1.24 - (13*k_u/12)), 0.85)
        
            #Final moment capacity##
            Ast_final = (self.des_bending_moment*1000*1000)/(phi*self.f_sy*lever_arm_final)
            # Ast_final = (self.des_bending_moment*1000*1000)/(0.85*500*lever_arm_assumed)

            return Ast_final
    
    def get_Icr(self):

        if (self.check_shear() == 'Increase size of beam, or increase ligatures/decrease spacing') or (self.ult_tensile_reo() == "Increase beam depth"):
            return "Beam Depth Not Adequate"
        else:
            n = self.steel_modulus/self.get_Ec()
            B = self.b/(n*self.ult_tensile_reo())
            d_tensile = 0.9*self.d
            d_comp = 72 ## Hard coded for now....##

            if self.As_top ==0:
                k =((2*d_tensile*B+1)**0.5-1)/(B*d_tensile)
                Icr = (self.b * (k*d_tensile)**3)/3 + n*self.ult_tensile_reo()*(self.d - k*self.d)**2
            if self.As_top !=0:
                r = ((n-1)*self.As_top)/(n*self.ult_tensile_reo())
                k = ((2*d_tensile*B*(1 + (r*d_comp/d_tensile)) + (1+r)**2)**0.5 - (1+r))/(B*d_tensile)

                Icr = (self.b * (k*d_tensile)**3)/3 + n*self.ult_tensile_reo()*(self.d - k*self.d)**2 + (n-1)*self.As_top*((k*d_tensile - d_comp)**2)

            return Icr
    

@dataclass
class conc_beam_service(concrete_rect_beam):
    time : Optional[int] = 10950

    def autogenous_shrinkage(self):
        
        if self.conc_grade <= 50.:
            autogenous_shrinkage = (0.07*self.conc_grade - 0.5)*50*(10**-6)*(1-math.exp(-0.07*self.time))
        elif self.conc_grade > 50.:
            autogenous_shrinkage = (0.08*self.conc_grade - 1.)*50*(10**-6)*(1-math.exp(-0.07*self.time))

        return autogenous_shrinkage

    def drying_shrinkage(self):

        k4 = 0.65
        basic_drying_shrinkage = (0.9-0.005*self.conc_grade)*800*10**(-6)

        th = (2*self.b*self.d)/(2*(self.b + self.d)) ## Have to check this calcualtion##
        alpha_1 = 0.8 + 1.2*math.exp(-0.005*th)
        k1 = (alpha_1*(self.time**0.8))/(self.time**0.8 + (0.15*th))

        drying_shrinkage = k1 * k4 * basic_drying_shrinkage

        return drying_shrinkage
    
    def final_shrinkage(self):
        return self.autogenous_shrinkage() + self.drying_shrinkage()
    
    def shrinkage_tensile_stress(self):

        reo_ratio = (self.ult_tensile_reo())/(self.b*self.d)

        shrinkage_tensile_stress = ((2.5*reo_ratio)/(1+50*reo_ratio))*200*1000*self.final_shrinkage()

        return shrinkage_tensile_stress

    def get_Mcr(self):
        section_modulus = (self.get_Ig())/(self.d/2)
        fctf = 0.6*(self.conc_grade)**0.5
        # shrinkage = self.shrinkage_tensile_stress()
        Mcr = section_modulus * (fctf)

        return Mcr/(1000000)
    
    def get_Ief_refined(self):
        Icr = self.get_Icr()
        Ig = self.get_Ig()
        Mcr = self.get_Mcr()
        Ms = self.ser_design_moment
        

        if self.ult_tensile_reo()=="Increase beam depth":
            return "Increase beam depth"
        else:
            reo_ratio = (self.ult_tensile_reo()+self.As_top)/(self.b*self.d)
            if reo_ratio<0.005:
                Ief = min(Icr/(1-((1-Icr/Ig)*(Mcr/Ms)**2)), 0.6*Ig)
            elif reo_ratio>= 0.005:
                Ief = min(Icr/(1-((1-Icr/Ig)*(Mcr/Ms)**2)), Ig)
            return Ief
    
    # def get_rigidity(self):
    #     Ief = self.get_Ief()
    #     Ec = self.get_Ec()
    #     rigidity_short_term = self.ser_design_moment/(Ec*Ief)
    #     ri

@dataclass
class conc_beam_shear(conc_beam_service):

    at_support : Optional[bool] = False
    lig_size: Optional[float] = 0.0
    lig_legs : Optional[float] = 0.0
    lig_spacing : Optional[float] = 0.0
    fsy_f : Optional[float] = 500.0
    des_shear : Optional[float] = 0.


    def get_Asv(self):

        area_lig_leg = (math.pi*(self.lig_size**2))/4
        Asv = (self.lig_legs*area_lig_leg)

        return Asv

    def get_Asvmin(self):

        Asvmin = max(0.06*(self.conc_grade)**0.5*self.b*self.lig_spacing/self.fsy_f , 0.35*self.b*self.lig_spacing/self.fsy_f)
        return Asvmin
    
    def get_Vuc(self):
        d_o = 0.9*self.d
        fcv = min(self.conc_grade**(1/3), 4.0)
        
        if self.get_Asv()>self.get_Asvmin():
            beta_1 = max(1.1*(1.6- (d_o/1000)), 1.1)
        else:
            beta_1 = max(1.1*(1.6- (d_o/1000)), 0.8)

        beta_2 = 1.0

        if self.at_support ==True:
            beta_3 = 2.0
        else:
            beta_3 = 1.0
        V_uc = beta_1*beta_2*beta_3*self.b*d_o*fcv*(((self.ult_tensile_reo())/(self.b*d_o))**(1/3))

        return V_uc/1000

    def get_Vus(self):
        d_o = 0.887*self.d
        Vus = (self.get_Asv()*self.fsy_f*d_o/self.lig_spacing)
        return Vus/1000
    
    def get_Vumin(self):
        d_o = 0.887*self.d
        Vu_min = max(self.get_Vuc() + (0.1*(self.conc_grade)**0.5*self.b*d_o)/1000, self.get_Vuc()+(0.6*self.b*d_o/1000))

        return Vu_min


    def get_Vumax(self):
        d_o = 0.887*self.d

        Vumax = 0.2*self.conc_grade*self.b*d_o/1000

        return Vumax
    
    def get_Vult(self):

        if self.ult_tensile_reo()=="Increase beam depth":
            return "Increase beam depth"
        else:
            V_u = self.get_Vuc() + self.get_Vus()
            V_umin = self.get_Vumin()
            V_umax = self.get_Vumax()

            if (V_u <= V_umax) and (V_u<=V_umin):
                return 0.7*V_umin   
            elif (V_u<=V_umax) and (V_u>=V_umin):
                return 0.7*V_u
    
        
    def check_shear(self):

        V_ult = self.get_Vult()
        if V_ult == "Increase beam depth":
            return "Increase beam depth"
        elif (self.des_shear > V_ult) or (self.des_shear>self.get_Vumax()):
            return "Increase size of beam, or increase ligatures/decrease spacing"
        elif self.des_shear <= V_ult:
            return f"The Section is working at {round((self.des_shear/V_ult)*100, 1)} %"
        
    def calc_shear_reo_weight(self):

        ##Calculate horizontal legnth of ligs##
        outer_lig_length = self.b - 60 - self.lig_size
        if self.lig_legs ==2 :
           lig_horizontal_length =  outer_lig_length
        elif self.lig_legs >2 and self.lig_legs % 2 == 0:
            s = outer_lig_length/(self.lig_legs -1)
            lig_horizontal_length = outer_lig_length
            inner_lig_length_horizontal = 0.
            for i in range(int(self.lig_legs/2 -1)):
                lig_horizontal_length = lig_horizontal_length - 2*s
                inner_lig_length_horizontal += lig_horizontal_length
            lig_horizontal_length = outer_lig_length + inner_lig_length_horizontal
        ## Calculate lig vertical length##
        lig_vertical_length = self.lig_legs*(self.d - 60)

        ## total lig weight per metre##
        weight_ligs = ((lig_horizontal_length+lig_vertical_length)/1000)*(math.pi*(self.lig_size**2)/(4*1000*1000))*(7850)*(1000/self.lig_spacing)
        return weight_ligs
      

def beam_deflections(span : float, concrete_beam, load_combo_list = []):
    if (concrete_beam.check_shear() == 'Increase size of beam, or increase ligatures/decrease spacing') or (concrete_beam.ult_tensile_reo() == 'Increase beam depth'):
        return "Beam Depth Not Adequate"
    else:
        ##Get Deflection##
        short_term_deflection = (5*load_combo_list[4]*(span*1000)**4)/(384*concrete_beam.get_Ief_refined()*concrete_beam.get_Ec()) + (load_combo_list[2]*1000*(span*1000)**3)/(48*concrete_beam.get_Ec()*concrete_beam.get_Ief_refined())
        sustained_deflection = (5*load_combo_list[5]*(span*1000)**4)/(384*concrete_beam.get_Ief_refined()*concrete_beam.get_Ec()) + (load_combo_list[3]*1000*(span*1000)**3)/(48*concrete_beam.get_Ec()*concrete_beam.get_Ief_refined())

        total_deflection = short_term_deflection + (concrete_beam.get_kcs())*sustained_deflection

        beam_analysis_dict = {}
        beam_analysis_dict["Short Term Deflection"] = round(short_term_deflection, 1)
        beam_analysis_dict["Total Deflection"] = round(total_deflection, 1)
        return beam_analysis_dict


def beam_deflections_refined(span : float, concrete_beam, load_combo_list = []):
    ##Bending moments##
    M_short_term_pl = (load_combo_list[2] * span *1000*1000)/4 
    M_short_term_udl = (load_combo_list[4]* (span*1000)**2)/8
    M_long_term_pl = (load_combo_list[3] * span *1000*1000)/4 
    M_long_term_udl = (load_combo_list[5]* (span*1000)**2)/8


    ## Deflection Short term##
    short_term_deflection_pl = (((span*1000)**2/48)*(4*M_short_term_pl))/(concrete_beam.get_Ec()*concrete_beam.get_Ief_refined())
    short_term_deflection_udl= (((span*1000)**2/96)*(10*M_short_term_udl ))/((concrete_beam.get_Ec()*concrete_beam.get_Ief_refined()))

    short_term_deflection = short_term_deflection_pl + short_term_deflection_udl

    ## Deflection Long Term##
    long_term_deflection_pl = (((span*1000)**2/48)*(4*M_long_term_pl))/(concrete_beam.get_Ec()*concrete_beam.get_Ief_refined())
    long_term_deflection_udl= (((span*1000)**2/96)*(10*M_long_term_udl))/((concrete_beam.get_Ec()*concrete_beam.get_Ief_refined()))

    sustained_deflection = long_term_deflection_pl + long_term_deflection_udl
    
    ##Get Deflection##
    total_deflection = short_term_deflection + (concrete_beam.get_kcs())*sustained_deflection

    beam_analysis_dict = {}
    beam_analysis_dict["Short Term Deflection"] = round(short_term_deflection, 1)
    beam_analysis_dict["Total Deflection"] = round(total_deflection, 1)
    beam_analysis_dict["Load Combo"] = load_combo_list

    return beam_analysis_dict


def calculate_reo_rate(concrete_beam):
    b = concrete_beam.b
    d = concrete_beam.d
    volume = b*d/(1000000)
    tenile_reo_weight = concrete_beam.ult_tensile_reo()/1000000 * 7850
    comp_reo_weight = (concrete_beam.As_top/(1000*1000))*7850

    shear_reo_weight = concrete_beam.calc_shear_reo_weight()

    return round((shear_reo_weight+tenile_reo_weight+comp_reo_weight)/volume , 1)

def reo_rate_comparison(concrete_beam):
    comparison_parameter = round(calculate_reo_rate(concrete_beam)/185.2 *100 ,1)

    if comparison_parameter == 100:
        return 0
    elif comparison_parameter > 100:
        return comparison_parameter
    elif comparison_parameter <100:
        return -1*comparison_parameter
    

    



    

