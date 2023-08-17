using Distributions, Plots

# Parameters
α, β = 5.0, 1.0  # α is the shape parameter, β is the rate parameter
n_K = 30 #dias que se usan hacia adelante que se usan para construir costos esperados
T = 360 #dias de simulacion
error_tol = 0.10 #error que tolera al limite, es mayor a 2% porque normalmente esta por encima por lo que la probabilidad baja.

lambdas = rand(Gamma(α, β), T + n_K + 1) #Las lambdas que caracterisan las distribuciones poisson de las demandas
M = Poisson.(lambdas)#Todas las distribuciones de demandas por dia

I_tilde = zeros(T + n_K + 1)#El numero tal que si el inventario es menor o igual que esto debemos rellenar
for t = 1:T + n_K + 1
    temp = round(Int, mean(M[t]))
    while true
        if 1 - cdf(M[t], temp) < error_tol
            break
        end
        temp += 1
    end
    I_tilde[t] = temp - 1
end

"Funcion que calcula el costo esperado si se hace restock el dia t y se coloca inventario I"
function loss_costo(I, t)
    #Obtiene distribucion de N (dias para detenerse)
    P_N = zeros(n_K)
    CDF_N = zeros(n_K)
    CDF_N[1] = 1 - cdf(M[t], I - I_tilde[t+1] - 1)
    P_N[1] = CDF_N[1]
    for k = 1:n_K-1
        #suma de poisson es poisson.
        #Esto es una aproximacion, porque la distribucion deberia tomar en cuenta la dependencia
        #En realidad se suman poissones truncados en I_tilde para t=0,... t-k-1. con una poisson completa para t=k. Pero Esto
        #es mas complejo de lo que quiero programar ahorita. Para resolverlo voy a normalizar la distribucion a 1, hay que tomar
        # en cuenta que estoy dandole peso a las n mas grandes.
        sum_distribution = Poisson(sum(lambdas[t:t+k]))
        CDF_N[1+k] = (1 - cdf(sum_distribution, (I - I_tilde[t+k+1] - 1)))
        CDF_N[1 + k] = max(CDF_N[1+k], CDF_N[k])
        P_N[1+k] = CDF_N[1+k] - CDF_N[k]
    end

    #P_N = P_N ./ sum(P_N)

    #First term
    E_1 = 0
    for N = 0:n_K-1
        E_1 += 100/(N+1) * P_N[N + 1]
    end

    #Second term
    E_2 = I
    

    #third term
    E_3 = 0
    for N = 0:n_K-1
        demand_sum_N = 0
        for k = 0:N
            demand_sum_N += mean(M[t + k]) * (N-k)/(N+1)
        end

        E_3 += P_N[N+1]*(demand_sum_N)
    end

    loss = E_1 + E_2 + E_3
    return loss

end

p1 = plot(collect(1:50), loss_costo.(collect(1:50), 1),
title="Costos diarios esperados al rellenar", xlabel="Inventario seleccionado", ylabel="\$", legend=false,
titlefont = 12, size = (500, 500))
savefig(p1, "plotcaso1.png")

"funcion que calcula el optimo I en tiempo t"
function restock_I!(t)
    rango_busqueda = collect(1:50) #rango donde busca minimo optimo.
    costos_I = loss_costo.(rango_busqueda, t)
    global I = argmin(costos_I)
    costos_acumulados[t] += 100 #records fixed expenses
end

"funcion que checa si se tiene que hacer restock"
function needs_restock(I, t)
    if 1-cdf(M[t], I) > error_tol
        return true
    else
        return false
    end
end 

#Inicia simulacion
global t_sim = 1
global I = 0
costos_acumulados = zeros(T+1)
I_record = zeros(T)
restock_days = [1]

#primer dia restock
#restock_I!(t_sim)
#t_sim += 1

#resto de los dias
while t_sim <= T
    if needs_restock(I, t_sim) #checks if it needs a restock
        restock_I!(t_sim)
        append!(restock_days, t_sim)
    end
    costos_acumulados[t_sim] += I #variable expenses are recorded
    m = rand(M[t_sim])
    global I = max(I - m, 0) #inventory is depleted
    I_record[t_sim] = I #records inventory level
    costos_acumulados[t_sim + 1] += costos_acumulados[t_sim]
    global t_sim += 1 #next day
    
end

histogram(restock_days[3:end] - restock_days[2:end-1])
plot(I_record)
sum(I_record .== 0)/T*100

#simulacion de estrategia no optima
global t_sim = 1
global I = 0
costos_acumulados2 = zeros(T+1)
I_record2 = zeros(T)
restock_days2 = [1]

while t_sim <= T
    if mod(t_sim, 7) == 1
        sum_distribution = Poisson(sum(lambdas[t_sim:t_sim+6]))
        global I = quantile(sum_distribution, 0.94)
        append!(restock_days, t_sim)
        costos_acumulados2[t_sim] += 100 #records fixed expenses
    end
    costos_acumulados2[t_sim] += I #variable expenses are recorded
    m = rand(M[t_sim])
    global I = max(I - m, 0) #inventory is depleted
    I_record2[t_sim] = I #records inventory level
    costos_acumulados2[t_sim + 1] += costos_acumulados2[t_sim]
    global t_sim += 1 #next day
end

sum(I_record2 .== 0)/T*100

p2 = plot(costos_acumulados, label = "estrategia: óptimo")
plot!(1:T, 100*collect(1:T), label = "estrategia: rellenar diario")
plot!(costos_acumulados2, label = "estrategia: rellenar cada semana", 
size = (500, 500),  xlabel="día", ylabel="\$", title = "costos acumulados")
savefig(p2, "plotcaso2.png")