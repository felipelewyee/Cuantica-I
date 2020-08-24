# Ejercicio 8. Hartree-Fock-Roothan

Tenemos que resolver:
\begin{equation}
\textbf{F} \textbf{C} = \textbf{S}\textbf{C}\varepsilon
\end{equation}

- <span style="color:green"> **Paso 1.** Especificar molécula: Coordenadas de los núcleos $\{R_A\}$, Carga de los núcleos $\{Z_A\}$, Número de electrones $(N)$, y funciones base $\{\phi_i\}$</span>
- <span style="color:green"> **Paso 2.** Calcular $S$, $H$, $(ij|kl)$.</span>
- <span style="color:green"> **Paso 3.** Proponer una matriz $C$.</span>
- <span style="color:black"> **Paso 4.** Calcular $P$, $J$ y $K$.</span>
- <span style="color:black"> **Paso 5.** Calcular $F=H+2J-K$</span>
- <span style="color:black"> **Paso 6.** Resolver $FC=SC\varepsilon$</span>
- <span style="color:black"> **Paso 7.** $E_{elec} = \sum_\mu \sum_\nu P_{\nu \mu} (H_{\mu \nu} + F_{\mu \nu})$<span>
- <span style="color:black"> **Paso 8.** ¿$E_i=E_{i-1}$?, Sí: acabé. No: volver a paso 4.</span>
- <span style="color:black">**Paso 9.** Calcular energía nuclear y sumarla a la energía electrónica.</span>



**Paso 1.** Especificar molécula: Coordenadas de los núcleos <span style="color:red">$\{R_A\}$</span>, Carga de los núcleos <span style="color:red">$\{Z_A\}$</span>, Número de electrones <span style="color:red">$(N)$</span>, y funciones base <span style="color:red">$\{\phi_i\}$</span>.

H2 = psi4.geometry("""




""")

psi4.set_options({'basis':''})

**Paso 2.** Calcular $S$, $H$, $(ij|kl)$.

<span style="color:red">
\begin{equation}
S_{ij} = \int \psi_i^*(r) \psi_j(r) dr
\end{equation}
</span>
\begin{equation}
H_{ij} = \int \psi_i^*(r) \hat{H} \psi_j(r) dr
\end{equation}
<span style="color:red">
\begin{equation}
H_{ij} = T_{ij} + V_{ij}
\end{equation}
</span>
\begin{equation}
\psi_i(r) = \sum_\mu a_\mu \phi_\mu(r)
\end{equation}
\begin{equation}
(\mu \nu|\sigma \lambda) = \int \int \phi_\mu^*(r_1) \phi_\nu(r_1) \frac{1}{r_{12}} \phi_\sigma^*(r_2) \phi_\lambda(r_2) dr_1 dr_2
\end{equation}

wfn = psi4.core.Wavefunction.build(H2)
mints = psi4.core.MintsHelper(wfn.basisset())

S = np.asarray(mints.ao_overlap())
print("----------------Matriz S----------------")
sp.pprint(sp.Matrix(S))

T = np.asarray(mints.ao_kinetic())
V = np.asarray(mints.ao_potential())
H=
print("----------------Matriz H----------------")
sp.pprint(sp.Matrix(H))

I = np.asarray(mints.ao_eri())
#print("Integrales (ij|kl):")
#print(I)

nbf = S.shape[0]
ndocc = wfn.nalpha()
#print(nbf)
#print(ndocc)

**Paso 3.** Proponer una matriz C.

C = 
#print(C)

- **Paso 4.** Calcular $P$, $J$ y $K$.
\begin{equation}
P_{\mu \nu} = \sum_a^{N/2} C_{\mu a} C_{\nu a}^*
\end{equation}
\begin{equation}
J_{\mu \nu} = \sum_{\lambda \sigma} P_{\lambda \sigma} (\mu \nu | \sigma \lambda)
\end{equation}
\begin{equation}
K_{\mu \nu} = \sum_{\lambda \sigma} P_{\lambda \sigma} (\mu \lambda | \sigma \nu)
\end{equation}
- **Paso 5.** Calcular $F=H+2J-K$
- **Paso 6.** Resolver $FC=SC\varepsilon$
- **Paso 7.** $E_{elec} = \sum_\mu \sum_\nu P_{\nu \mu}(H_{\mu \nu} + F_{\mu \nu})$
- **Paso 8.** ¿$E_i=E_{i-1}$?, Sí: acabé. No: volver a paso 4.



while():
    print("---------------------------------ITERACION---------------------------------")  

    print("------------Matriz C Entrada------------")
    sp.pprint(sp.Matrix(C))    
    
    #Paso 4. Calcular P, J y K
    P = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(ndocc):
                    P[i][j] = P[i][j] + C[i][k]*C[j][k]
    print("----------------Matriz P----------------")
    sp.pprint(sp.Matrix(P))
                    
    J = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(nbf):
                for l in range(nbf):
                    J[i][j] = J[i][j] + P[k][l]*I[i][j][l][k]
    print("----------------Matriz J----------------")
    sp.pprint(sp.Matrix(J))
                    
    K = np.zeros((nbf,nbf))
    for i in range(nbf):
        for j in range(nbf):
            for k in range(nbf):
                for l in range(nbf):
                    K[i][j] = K[i][j] + P[k][l]*I[i][l][k][j]         
    print("----------------Matriz K----------------")
    sp.pprint(sp.Matrix(K))                    

    #Paso 5. Calcular F = H + 2J - K
    F = 
    print("----------------Matriz F----------------")
    sp.pprint(sp.Matrix(F))
    
    #Paso 6. Resolver FC=SCE
    E,C = 
    print("-------------Matriz C Salida-------------")
    sp.pprint(sp.Matrix(C))
    
    #Paso 7. Calcular E=sum_i sum_j P_ji (H_ij+F_ij)
    E_elec = 0.0
    for i in range(nbf):
        for j in range(nbf):
            E_elec = E_elec + P[j][i]*(H[i][j] + F[i][j])

    
    print("Energia Electronica: ",E_elec)

    #Paso 8. ¿$E_i=E_{i-1}$?, Si: acabe. No: volver a paso 4.

        

**Paso 9.** Calcular energía nuclear y sumarla a la energía electrónica.
\begin{equation}
E_{\rm nuc} = \sum_{A>B}\frac{Z_AZ_B}{r_{AB}}
\end{equation}
\begin{equation}
E_{Tot} = E_{elec} + E_{nuc}
\end{equation}

E_nuc = H2.nuclear_repulsion_energy() 
print("Energia nuclear: ", E_nuc)
E_T = 
print("Energia Total: ", E_T)

# Método Simple. PSI4

Geometría del agua

|Átomo|X (A)|Y (A)|Z (A)|
|----|----|----|----|
|H|0.0000|0.0000|0.0000|
|H|0.0000|0.0000|0.7414|








# Referencias

- Szabo, A.; Ostlund, N. S. Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory; Dover Publications: Mineola, N.Y, 1996.

Hecho por Juan Felipe Huan Lew Yee y Jorge Martín del Campo Ramírez para la clase de Química Cuántica I del ciclo 2019-I. Universidad Nacional Autónoma de México.

Este archivo puede distribuise libremente y ser considerado Open Source. Si deseas modificarlo para su distribución, solo se pide conservar el nombre de los autores originales.