import streamlit as st
# Запуск: streamlit run d:\sml\sml\sml.py
# Выход: Ctrl+C

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import *
import altair as alt

st.title('Построение эллипса фон Мизеса, расчет параметров и характеристик обсадной трубы')

# Исходные данные:
#D = 177.8 # Диаметр трубы, мм
D = st.number_input('Диаметр трубы, мм', value=177.8, step=None)
#t = 9.19 # Толщина стенки, мм
t = st.number_input('Толщина стенки, мм', value=9.19, step=None)

data_pipe = { 'Диаметр обсадной трубы, мм' : ['101,60', '114,30','127,00', '139,70', '146,05', '168,28', '177,80', '193,68', '219,08', '244,48', '250,83', '273,05', '298,45', '323,85', '339,72', '350,52', '376,76', '406,40', '425,45', '473,08', '508,00'], 'Варианты толщин стенок, мм': ['6,50', '5,21 / 5,69 / 6,35 / 7,37 / 8,56 / 10,20', '5,59 / 6,43 / 7,52 / 9,19 / 10,70 / 11,10 / 12,14 / 12,70', '6,20 / 6,98 / 7,72 / 9,17 / 10,54 / 12,70', '6,50 / 7,00 / 7,70 / 8,50 / 9,50 /10,70', '7,32 / 8,94 / 10,59 / 12,06', '5,87 / 6,91 / 8,05 / 9,19 / 10,36 / 11,51 /12,65 / 13,72 / 15,00', '7,62 / 8,33 / 9,52 / 10,92 / 12,70 / 14,27 / 15,11 / 15,88', '6,17 / 7,72 / 8,94 / 10,16 / 11,43 / 12,70 / 14,15', '7,92 / 8,94 / 10,03 / 11,05 / 11,99 /13,84 / 15,11', '15,88', '7,09 / 8,89 / 10,16 / 11,43 / 12,57 / 13,84 / 15,11 / 16,50', '8,46 / 9,53 / 11,05 / 12,42 / 13,56 / 14,78', '7,70 / 8,50 / 9,50 / 11,00 / 12,40 / 14,00', '8,38 / 9,65 / 10,92 / 12,19 / 13,06 / 14,00 /15,40', '9,00 / 10,00 / 11,00 / 12,00', '9,00 / 10,00 / 11,00 / 12,00', '9,53 / 11,13 / 12,57 / 16,66', '10,00 / 11,00', '11,05', '11,13 / 12,70 / 16,13']}
frame_pipe = pd.DataFrame(data_pipe)
st.table(frame_pipe)

#Sy = 379 # Предел текучести, МПа
Sy = st.number_input('Минимальный предел текучести, МПа', value=379, step=None)

data_Sy = {'Группа прочности' : ['H40', 'J55 (K55, Д)', 'K72 (К)', 'N80 (N80-Q, L80, Е)', 'C90', 'C95 (T95, Л)','P110 (М)', 'Q125', 'Q135 (Р)', '- 140', '- 150' ], 'Минимальный предел текучести, МПа' :[276, 379, 490, 552, 621, 655, 758, 862, 930, 966, 1035]}
frame_Sy = pd.DataFrame(data_Sy)
#st.write(frame_Sy)
st.table(frame_Sy)

#K_eff_min = 50 # Эффективность соединения на сжатие, %
K_eff_min = st.number_input('Эффективность соединения на сжатие, %', value=50, step=None)
#K_eff_pls = 100 # Эффективность соединения на растяжение, %
K_eff_pls = st.number_input('Эффективность соединения на растяжение, %', value=100, step=None)

# Коэффициенты запаса для построения элипса фон Мизеса с учетом коэфф.запаса
# Нужно доработать ввод коэффициентов запаса

d = D-2*t # Внутренний диаметр трубы, мм
st.write('Внутренний диаметр трубы, мм :', round(d, 2))
Ap = (math.pi/4)*(D**2-d**2) # Площадь поперечного сечения трубы, мм^2
st.write('Площадь поперечного сечения трубы, мм^2 :',  math.floor(Ap))

F_max = Sy*Ap # Максимальная растягивающая нагрузка до предела текучести, Н
st.write('Максимальная растягивающая нагрузка до предела текучести, кН :', math.floor(F_max/1000))

F_max_con = F_max*K_eff_min/100 # Максимальная нагрузка на сжатие для соединения с учетом эффективности соединения на сжатие, Н
st.write('Максимальная нагрузка на сжатие для соединения с учетом эффективности соединения на сжатие, кН :', math.floor(F_max_con/1000))

# Минимальное внутреннее давление до предела текучести, МПа
Pi_max = 2*Sy*0.875*t/D
#Pi_max = 2*Sy*1.000*t/D # Верхняя граница элипса. Коэффициент 1.000, вместо 0,875

st.write('Минимальное внутреннее давление до предела текучести, МПа :', round(Pi_max, 2))

# Функция расчета наружнего давления
def Collapse (x,dn,t,Ap,Sy) :
        # x-расчетная нагрузка, Н
        
        Syr = Sy/0.00689475729317831
        Sa = x/Ap/0.00689475729317831
        T = ((1-0.75*(Sa/Syr)**2)**0.5-0.5*Sa/Syr)*Syr
        
        c0=2.8762
        c1=0.10679/100000
        c2=0.21302/10000000000
        c3=-0.53132/10000000000000000
        c4=0.026233
        c5=0.50609/1000000
        c6=-465.93
        c7=0.030867
        c8=-0.10483/10000000
        c9=0.36989/10000000000000
        c10=46.95*1000000
        
        F1=c0+c1*T+c2*(T**2)+c3*(T**3)
        F2=c4+c5*T
        F3=c6+c7*T+c8*(T**2)+c9*(T**3)
        RF=F2/F1
        F4=c10*((3*RF/(2+RF))**3)/T/(3*RF/(2+RF)-RF)/((1-3*RF/(2+RF))**2)
        F5=F4*RF
        
        Dtyp=(((((F1-2)**2)+8*(F2+F3/T))**0.5)+(F1-2))/(2*(F2+(F3/T))) #t1
        Dtte=(2+(F2/F1))/((3*F2)/F1) #t3
        Dtp=(T*(F1-F4))/(F3+T*(F2-F5)) #t2
        
        P_yp = 2*T*((dn/t-1)/((dn/t)**2))
        P_p = T*((F1/(dn/t)-F2))-F3
        P_t = T*(F4/(dn/t)-F5)
        P_e = (46.95*(10)**6)/((dn/t)*((dn/t-1)**2))
        
        # Доработать расчет наружных давлений
        if dn/t < Dtyp:
            Pcr=P_yp
        else:
            if np.any((Dtyp < dn/t) & (dn/t < Dtp)):
                Pcr=P_p
            else:
                if np.any((Dtp < dn/t) & (dn/t < Dtte)):
                    Pcr=P_t
                else:
                    Pcr=P_e
        return Pcr*0.00689475729317831 #Значение наружнего давления с отрицательным знаком в нижних квадрантах

st.write('Минимальное наружнее давление до предела текучести по стандарту ISO 10400, МПа :', round(Collapse (0,D,t,Ap,Sy), 2))

q = 1.017*7800*(math.pi/4*(D**2-d**2))/10**8

st.write('Вес одного метра колонны (коэффициент к массе трубы на муфту 1,017), кН :', round(q,3))

Q_k = st.number_input('Длина колонны, м', value=1000, step=None)

st.write('Вес колонны в водухе, кН :', round(q*Q_k,2))
   
def Pi1 (F,D,t,Ap,Sy):
        d = D-2*t
        Fa = F
        Sa = Fa/Ap
        k_pi = (D**2+d**2)/(D**2-d**2)
        k_c = Sa**2-Sy**2
        k_b = (1-k_pi)*Sa
        k_a = k_pi**2+k_pi+1
        Pi = (-k_b+np.sqrt(k_b**2-4*k_a*k_c))/(2*k_a) # ...-k_b+np.sqrt... тут знак "+"
        return Pi

#st.write('Минимальное наружнее давление до предела текучести по стандарту ISO 10400, МПа :', round(Pi1(0,D,t,Ap,Sy), 2))


#  Критерий текучести тела трубы по фон Мизесу при действии внутреннего и наружного давлений и осевого напряжения по стандарту 
# Функция создания массива данных для построения графика 1 и 2 квадрантов. Нагрузка/внутреннее давление
def Pi (D,t,Sy):
    Pi_F = np.array([[],[]])
    d = D-2*t
    Ap = (math.pi/4)*(D**2-d**2)
    F_max = Sy*Ap
       
    # Коэффициент k_pi
    k_pi = (D**2+d**2)/(D**2-d**2)
        
    # Нагрузка в точке смены знака +/- при определении внутреннего давления
    F_max_kr = (2*Sy*Ap*np.sqrt(k_pi**2+k_pi+1))/(np.sqrt(3)*(k_pi+1))
    F_max_kr = math.trunc(F_max_kr) # Отбросим дробную часть, что-бы избежать вызова исключений при вычислении квадратного корня
    
    # Массив значений нагрузки от -F_max до F_max_kr. 1000 точек
    F_arr1 = np.linspace(-F_max,F_max_kr,1000)
    # Массив значений нагрузки от F_max_kr до F_max. 200 точек
    F_arr2 = np.linspace(F_max_kr,F_max,200)
    
    # Функция вычисления внутренних давлений квадрант 1  в интервале нагрузок от -F_max до F_max_kr
    def Pi1 (F,D,t,Ap,Sy):
        d = D-2*t
        Fa = F
        Sa = Fa/Ap
        k_pi = (D**2+d**2)/(D**2-d**2)
        k_c = Sa**2-Sy**2
        k_b = (1-k_pi)*Sa
        k_a = k_pi**2+k_pi+1
        Pi = (-k_b+np.sqrt(k_b**2-4*k_a*k_c))/(2*k_a) # ...-k_b+np.sqrt... тут знак "+"
        return Pi
    
    # Функция вычисления внутренних давлений квадрант 2 до в интервале нагрузок от F_max_kr до F_max
    def Pi2 (F,D,t,Ap,Sy):
        d = D-2*t
        Fa = F
        Sa = Fa/Ap
        k_pi = (D**2+d**2)/(D**2-d**2)
        k_c = Sa**2-Sy**2
        k_b = (1-k_pi)*Sa
        k_a = k_pi**2+k_pi+1
        Pi = (-k_b-np.sqrt(k_b**2-4*k_a*k_c))/(2*k_a) # ...-k_b-np.sqrt... тут знак "-"
        return Pi
    
    Pi_arr1 = np.array([Pi1(F,D,t,Ap,Sy) for F in F_arr1]) # Расчет массива давлений Pi1
    
    Pi_arr2 = np.array([Pi2(F,D,t,Ap,Sy) for F in F_arr2]) # Расчет массива давлений Pi2
    
    Pi_F = [np.concatenate((F_arr1, F_arr2)),np.concatenate((Pi_arr1, Pi_arr2))] # Объединение массивов нагрузок и давлений
          
    return Pi_F # Вывод

Pi_plot = Pi(D,t,Sy)

# Функция создания массива данных для построения графика 1 и 2 квадрантов. Нагрузка/наружнее давление
def Po (D,t,Sy):
    Po_F = np.array([[],[]])
    d = D-2*t
    Ap = (math.pi/4)*(D**2-d**2)
    F_max = Sy*Ap
    
    # Коэффициент k_po
    k_po = (2*D**2)/(D**2-d**2)
        
    # Нагрузка в точке смены знака +/- при определении внутреннего давления
    F_max_kr = 2*Sy*Ap/np.sqrt(3)
    F_max_kr = math.trunc(F_max_kr) # Отбросим дробную часть, что-бы избежать вызова исключений при вычислении квадратного корня
    
    # Массив значений нагрузки от -F_max до F_max_kr. 200 точек
    F_arr1 = np.linspace(-F_max,-F_max_kr,200)
    # Массив значений нагрузки от F_max_kr до F_max. 1000 точек
    F_arr2 = np.linspace(-F_max_kr,F_max,1000)
    
    # Функция вычисления наружнего давлений квадрант 3 в интервале нагрузок от -F_max_kr до F_max
    def Po3(F,D,t,Ap,Sy):
        d = D-2*t
        Fa = F
        Sa = Fa/Ap
        k_po = (2*D**2)/(D**2-d**2)
        k_c = Sa**2-Sy**2
        k_b = k_po*Sa
        k_a = k_po**2
        Po = (-k_b-(k_b**2-4*k_a*k_c)**0.5)/(2*k_a) # ...-k_b-(k_b**2... тут знак "-"
        return -Po #Значение наружнего давления с отрицательным знаком в нижних квадрантах
    
    # Функция вычисления наружнего давлений квадрант 3 в интервале нагрузок от -F_max_kr до F_max
    def Po4(F,D,t,Ap,Sy):
        d = D-2*t
        Fa = F
        Sa = Fa/Ap
        k_po = (2*D**2)/(D**2-d**2)
        k_c = Sa**2-Sy**2
        k_b = k_po*Sa
        k_a = k_po**2
        Po = (-k_b+(k_b**2-4*k_a*k_c)**0.5)/(2*k_a) # ...-k_b+(k_b**2... тут знак "+"
        return -Po #Значение наружнего давления с отрицательным знаком в нижних квадрантах
    
    Po_arr1 = np.array([Po3(F,D,t,Ap,Sy) for F in F_arr1]) # Расчет массива давлений Pi1
    
    Po_arr2 = np.array([Po4(F,D,t,Ap,Sy) for F in F_arr2]) # Расчет массива давлений Pi2
    
    Po_F = [np.concatenate((F_arr1, F_arr2)),np.concatenate((Po_arr1, Po_arr2))] # Объединение массивов нагрузок и давлений
          
    return Po_F # Вывод

Po_plot=Po(D,t,Sy)

#Расчет максимального наружного/сминающего давления в МПа по стандарту ISO 10400

def Po_ISO (D,t,Sy):
    
    d = D-2*t
    Ap = (math.pi/4)*(D**2-d**2)
    F_max = Sy*Ap
    
    # Функция расчета наружнего давления
    def Collapse (x,dn,t,Ap,Sy) :
        # x-расчетная нагрузка, Н
        
        Syr = Sy/0.00689475729317831
        Sa = x/Ap/0.00689475729317831
        T = ((1-0.75*(Sa/Syr)**2)**0.5-0.5*Sa/Syr)*Syr
        
        c0=2.8762
        c1=0.10679/100000
        c2=0.21302/10000000000
        c3=-0.53132/10000000000000000
        c4=0.026233
        c5=0.50609/1000000
        c6=-465.93
        c7=0.030867
        c8=-0.10483/10000000
        c9=0.36989/10000000000000
        c10=46.95*1000000
        
        F1=c0+c1*T+c2*(T**2)+c3*(T**3)
        F2=c4+c5*T
        F3=c6+c7*T+c8*(T**2)+c9*(T**3)
        RF=F2/F1
        F4=c10*((3*RF/(2+RF))**3)/T/(3*RF/(2+RF)-RF)/((1-3*RF/(2+RF))**2)
        F5=F4*RF
        
        Dtyp=(((((F1-2)**2)+8*(F2+F3/T))**0.5)+(F1-2))/(2*(F2+(F3/T))) #t1
        Dtte=(2+(F2/F1))/((3*F2)/F1) #t3
        Dtp=(T*(F1-F4))/(F3+T*(F2-F5)) #t2
        
        P_yp = 2*T*((dn/t-1)/((dn/t)**2))
        P_p = T*((F1/(dn/t)-F2))-F3
        P_t = T*(F4/(dn/t)-F5)
        P_e = (46.95*(10)**6)/((dn/t)*((dn/t-1)**2))
        
        # Доработать расчет наружных давлений
        if dn/t < Dtyp:
            Pcr=P_yp
        else:
            if np.any((Dtyp < dn/t) & (dn/t < Dtp)):
                Pcr=P_p
            else:
                if np.any((Dtp < dn/t) & (dn/t < Dtte)):
                    Pcr=P_t
                else:
                    Pcr=P_e
        return -Pcr*0.00689475729317831 #Значение наружнего давления с отрицательным знаком в нижних квадрантах
    
    F_max_kr = 2*Sy*Ap/np.sqrt(3)
    F_max_kr = math.trunc(F_max_kr) # Отбросим дробную часть, что-бы избежать вызова исключений при вычислении квадратного корня
    
    # Определение максимального значения функции в квадранте 4. Максимальная нагрузка для расчета
    F_kr_kr=minimize_scalar(lambda x: -Collapse(x,D/25.4,t/25.4,Ap,Sy), bounds=[-F_max_kr,F_max], method='bounded')
        
    # Массив значений нагрузки от F_max_kr до F_kr_kr.x. 1000 точек
    F_arr = np.linspace(-F_max_kr,F_kr_kr.x,1000)
    
    
    Po_ISO_arr = np.array([Collapse(x,D/25.4,t/25.4,Ap,Sy) for x in F_arr]) # Расчет массива давлений
    
    Po_ISO_F = [F_arr,Po_ISO_arr] # Объединение массивов нагрузок и давлений
    
    return  Po_ISO_F # Вывод

Po_ISO_plot = Po_ISO (D,t,Sy)

# Определение граничных значений для графика ограничения по сжатию соединения
Con_Pi = Pi1(-F_max_con,D,t,Ap,Sy)
Con_Po = Collapse(-F_max_con,D/25.4,t/25.4,Ap,Sy)

# Построение графика
plt.rc('figure', figsize=(15, 10)) # Размер графика
fig, ax = plt.subplots()
ax.grid(True, which='both')
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
plt.xlabel("Осевая нагрузка, кН",labelpad=270)
plt.ylabel("Давление, МПа",labelpad=420)
ax.yaxis.set_minor_locator(plt.MultipleLocator(10)) # Цена деления
ax.xaxis.set_minor_locator(plt.MultipleLocator(500)) # Цена деления

plt.plot(Pi_plot[0]/1000,Pi_plot[1],color='blue',label='Эллипс фон Мизеса') # Элипс нагрузка/внутреннее давление квадранты 1, 2
plt.plot(Po_plot[0]/1000,Po_plot[1],color='blue') # Элипс нагрузка/наружнее давление квадранты 3, 4
plt.plot(Po_ISO_plot[0]/1000,Po_ISO_plot[1],color='red',label='Наружнее давление определенное по ISO 10400') # Элипс нагрузка/наружнее давление квадранты 3, 4 по ISO 10400
plt.plot([-F_max_con/1000,-F_max_con/1000],[Con_Pi,-Con_Po],color='green',label='Ограничение по сжатию для соединения') #Ограничение по сжатию для соединения

# Элипс фон Мизеса с учетом коэффициентов запаса
#plt.plot(Pi_plot[0]/1000/1.25,Pi_plot[1]/1.1,color='black',label='Элипс фон Мизеса с учетом коэффициентов запаса')
#plt.plot(Po_plot[0]/1000/1.25,Po_plot[1]/1.125,color='black')
#plt.plot(Po_ISO_plot[0]/1000/1.25,Po_ISO_plot[1]/1.125,color='black')

#plt.fill_between(Pi_plot[0]/1000/1.25,Pi_plot[1]/1.1,color="b", alpha=0.1)
#plt.fill_between(Po_plot[0]/1000/1.25,Po_plot[1]/1.125,color="b", alpha=0.1)
plt.legend(loc='upper left') # Настройка расположения легенды

# Вывод графика
st.pyplot(plt)

st.header('Расчет изменения наружнего/внутреннего давления с изменением толщины стенки')
#PN  Давление граничное, МПа
PN = st.number_input('Граничное наружнее давление, МПа (по умолчанию 100% от максимальнрого наружнего давления)', value=1.00*Collapse (0,D,t,Ap,Sy), max_value=Collapse (0,D,t,Ap,Sy), step=None)
#st.number_input('Граничное давление, МПа', value=round((0,75*Collapse (0,D,t,Ap,Sy)), 2), step=None)
PV = st.number_input('Граничное внутреннее давление, МПа (по умолчанию 85% от максимальнрого внутреннего давления)', value=0.85*2*Sy*0.875*t/D, max_value=2*Sy*0.875*t/D, step=None)

sigma_t = st.number_input('Процент износа толщины стенки, % (по умолчанию 25% от номинальной толщины стенки, мксимальное значение 50%)', value=25, min_value=0, max_value=50, step=None)

st.header('Наружнее давление:')
t_ar = np.append( np.arange(-t, -t*(100-sigma_t)/100, 0.01), -t)*-1
P_n = np.array([Collapse(0,D/25.4,t/25.4,Ap,Sy) for t in t_ar])
#P_ras = np.full(t_ar.size, round(Collapse (0,D,t,Ap,Sy), 2))
P_ras = np.full(t_ar.size, PN)
source1 = pd.DataFrame({'t' : t_ar, 'P': P_n, 'P_ras': P_ras})

#alt_g = alt.Chart(source).mark_line(color='red').encode( x='t', y='P', tooltip=['t', 'P'])
#alt_g2 = alt.Chart(source).mark_line(color='blue').encode(y='P_ras')
alt_g = alt.Chart(source1).encode( x='t', tooltip=['t', 'P', 'P_ras'])
alt_g1 = alt_g.mark_line(color='red').encode(y='P')
#alt_g1 = alt_g.mark_line(color='red').encode(alt.Y('P', title='Давление граничное, МПа'))
alt_g2 = alt_g.mark_line(color='blue').encode(y='P_ras')
alt_g.properties(width=600)
st.altair_chart(alt_g1 + alt_g2, use_container_width=True)

st.header('Внутреннее давление:')
P_v = np.array([2*Sy*0.875*t/D for t in t_ar])
P_ras_v = np.full(t_ar.size, PV)

source2 = pd.DataFrame({'t' : t_ar, 'P': P_v, 'P_ras': P_ras_v})

alt_g_v = alt.Chart(source2).encode( x='t', tooltip=['t', 'P', 'P_ras'])
alt_g1_v = alt_g_v.mark_line(color='red').encode(y='P')
#alt_g1 = alt_g.mark_line(color='red').encode(alt.Y('P', title='Давление граничное, МПа'))
alt_g2_v = alt_g_v.mark_line(color='blue').encode(y='P_ras')
alt_g_v.properties(width=600)
st.altair_chart(alt_g1_v + alt_g2_v, use_container_width=True)