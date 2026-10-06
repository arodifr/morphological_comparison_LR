library(geomorph)
library(dplyr)
library(ggplot2)
library(patchwork)


#### TRAITS ####
# Escenario con información completa
TRAITS_FULL <- c(
  "outline",
  "cejas_forma",
  "cejas_tail",
  "cejas_size",
  "cejas_separacion",
  "ojos_forma",
  "ojos_inclinacion",
  "ojos_size",
  "ojos_separacion",
  "nariz_altura",
  "nariz_anchura",
  "philtrum_altura",
  "boca_forma",
  "boca_comisuras",
  "labios_size",
  "boca_anchura",
  "barbilla_altura"
)

# Escenarios con información reducida
TRAITS_CAP <- c(
  "outline",
  "nariz_anchura",
  "philtrum_altura",
  "boca_forma",
  "boca_comisuras",
  "labios_size",
  "boca_anchura",
  "barbilla_altura"
)

TRAITS_GLASSES <- c(
  "outline",
  "cejas_forma",
  "cejas_tail",
  "cejas_size",
  "cejas_separacion",
  "nariz_anchura",
  "philtrum_altura",
  "boca_forma",
  "boca_comisuras",
  "labios_size",
  "boca_anchura",
  "barbilla_altura"
)

TRAITS_MASK <- c(
  "cejas_forma",
  "cejas_tail",
  "cejas_size",
  "cejas_separacion",
  "ojos_forma",
  "ojos_inclinacion",
  "ojos_size",
  "ojos_separacion"
)

# Leyendo los datos
#Apuntar este directorio a donde están los archivos
slids <- read.table(file.path("./csliders_aberdeenIA.txt"), header = FALSE )
coord_raw <- read.table(file.path("./coords_aberdeenIA.txt"), header = FALSE )
clasi <- read.table( file.path("./class_aberdeenIA.txt"), header = TRUE )

colnames(coord_raw)[1] <- "ID"
datos <- clasi %>% left_join(coord_raw, by = "ID") %>% filter(sex == "M")

clasiS <- datos[, 1:4]
coordS <- datos[, 5:124]


## RASGOS FACIALES shape y size ####
n_lm <- 60
coord <- arrayspecs(coordS, p = n_lm, k = 2 )
face.gpa <- gpagen(coord, curves = slids, PrinAxes = FALSE, ProcD = TRUE )
coord.Procrus <- two.d.array(face.gpa$coords)
rownames(coord.Procrus) <- clasiS$name

## Funciones 
# Para seleccionar los landmarks de cada rasgo
get_lm <- function(coords, landmarks) {
  cols <- as.vector(
    rbind(
      2 * landmarks - 1,
      2 * landmarks
    )
  )
  coords[, cols, drop = FALSE]
}
# Para graficar las coordenadas y checar al selección
plot_GMcoords <- function(coords, dimension){
  num_lms <- dim(coords)[2]/2
  plotAllSpecimens(arrayspecs(coords, num_lms, dimension)) 
}
# Para obtener los PC scores
get_pca_scores <- function(coords, landmarks) {
  
  x <- get_lm(coords, landmarks)
  
  gpa <- gpagen(
    arrayspecs(x, ncol(x) / 2, 2)
  )
  
  pca <- gm.prcomp(gpa$coords)
  
  list(
    gpa = gpa,
    pca = pca,
    scores = as.data.frame(pca$x)
  )
}

## outline
plot_GMcoords(get_lm(coord.Procrus, c(2, 42:59)), 2)
outline_scores <- get_pca_scores(coord.Procrus, c(2, 42:59))
# shape
DFoutline <- data.frame(
  clasiS[, 1:4],
  cara = outline_scores$scores$Comp1
)

## eyes
plot_GMcoords(get_lm(coord.Procrus, 18:25), 2)
eye_scores <- get_pca_scores(coord.Procrus, 18:25)

DFeye <- data.frame(
  clasiS[, 1:4],
  forma = eye_scores$scores$Comp2,
  inclinacion = eye_scores$scores$Comp1
)

plot_GMcoords(get_lm(coord.Procrus, 22:25), 2)
eyeR_scores <- get_pca_scores(coord.Procrus, 18:21)
eyeL_scores <- get_pca_scores(coord.Procrus, 22:25)
DFeye$size <- (eyeR_scores$gpa$Csize + eyeL_scores$gpa$Csize)/2

## boca
plot_GMcoords(get_lm(coord.Procrus, c(3:12, 60)), 2)
mouth_scores <- get_pca_scores(coord.Procrus, c(3:12, 60))

DFmouth <- data.frame(
  clasiS[, 1:4],
  forma = mouth_scores$scores$Comp1,
  inclinacion = mouth_scores$scores$Comp2,
  size = mouth_scores$gpa$Csize
)

## cejas
plot_GMcoords(get_lm(coord.Procrus, c(26:27, 36:41)), 2)
eyebrowR_scores <- get_pca_scores(coord.Procrus, c(26:27, 36:41))
eyebrowL_scores <- get_pca_scores(coord.Procrus, c(29:28, 30:35))

DFeyebrows <- data.frame(
  clasiS[, 1:4],
  forma = eyebrowR_scores$scores$Comp1, # rectangular gruesa o delgada
  tail = eyebrowR_scores$scores$Comp2, # cola en punta o cuadrada
  size = (eyebrowR_scores$gpa$Csize + eyebrowL_scores$gpa$Csize)/2
)


# RASGOS FACIALES distancias ####
# Función
euclidean <- function(a, b) sqrt(sum((a - b)^2))


endo_endo <- c()
cejas <- c()
glabella <-list()
glabela_subnasal <-c()
alar_alar <- c()
chelion_chelion <-c()
subnasal_labSup <- c()
nation_labInf <-c()
glabela_nation <- c()

for (i in 1:dim(coord.Procrus)[1]) {
  glabella[[i]] <- ((coord.Procrus[i,(27*2-1):(27*2)] + coord.Procrus[i,(29*2-1):(29*2)]) / 2)
}

for (i in 1:dim(coord.Procrus)[1]) {
  endo_endo[i] <-euclidean(coord.Procrus[i,(18*2-1):(18*2)], coord.Procrus[i,(22*2-1):(22*2)])
  cejas[i] <-euclidean(coord.Procrus[i,(27*2-1):(27*2)], coord.Procrus[i,(29*2-1):(29*2)])
  glabela_subnasal[i] <-euclidean(as.numeric(glabella[[i]]), coord.Procrus[i,(13*2-1):(13*2)])
  alar_alar[i] <-euclidean(coord.Procrus[i,(14*2-1):(14*2)], coord.Procrus[i,(15*2-1):(15*2)])
  chelion_chelion[i] <-euclidean(coord.Procrus[i,(7*2-1):(7*2)], coord.Procrus[i,(6*2-1):(6*2)])
  subnasal_labSup[i] <-euclidean(coord.Procrus[i,(8*2-1):(8*2)], coord.Procrus[i,(13*2-1):(13*2)])
  nation_labInf[i] <-euclidean(coord.Procrus[i,(2*2-1):(2*2)], coord.Procrus[i,(3*2-1):(3*2)])
  glabela_nation[i] <-euclidean(glabella[[i]], coord.Procrus[i,(2*2-1):(2*2)])
}


DFojosSep <- data.frame(clasiS[,1:4], 
                        ojosSep = endo_endo,
                        cejasSep = cejas)
DFnariz <- data.frame(clasiS[,1:4], 
                      altura = glabela_subnasal,
                      anchura = alar_alar) 
DFboca <- data.frame(clasiS[,1:4], anchura = chelion_chelion)
DFphiltrum <- data.frame(clasiS[,1:4], altura = subnasal_labSup)
DFbarbilla <- data.frame(clasiS[,1:4], altura = nation_labInf)


# DF DE TODOS LOS RASGOS FACIALES ####
# Shape y size
DFoutline
DFeyebrows
DFeye
DFmouth
# Distancias
DFojosSep # separacion cejas y ojos
DFnariz # altura y anchura
DFboca # anchura
DFphiltrum # altura
DFbarbilla # altura 

# DF
features <- clasiS[, 1:4]
features$outline <- DFoutline$cara
features$cejas_forma <- DFeyebrows$forma
features$cejas_tail <- DFeyebrows$tail
features$cejas_size <- DFeyebrows$size
features$cejas_separacion <- DFojosSep$cejasSep
features$ojos_forma <- DFeye$forma
features$ojos_inclinacion <- DFeye$inclinacion
features$ojos_size <- DFeye$size
features$ojos_separacion <- DFojosSep$ojosSep
features$nariz_altura <- DFnariz$altura
features$nariz_anchura <- DFnariz$anchura
features$philtrum_altura <- DFphiltrum$altura
features$boca_forma <- DFmouth$forma
features$boca_comisuras <- DFmouth$inclinacion
features$labios_size <- DFmouth$size
features$boca_anchura <- DFboca$anchura
features$barbilla_altura <- DFbarbilla$altura
head(features)


# ASIGNAR CATEGORIAS ####
## Funciones 
# Para obtener la media, SD del DF de referencia
get_cutpoints <- function(x, cutoff = 1) 
  {
  
  mu <- mean(x, na.rm = TRUE)
  s  <- sd(x, na.rm = TRUE)
  
  tibble(
    mean = mu,
    sd = s,
    lower = mu - (cutoff * s),
    upper = mu + (cutoff * s)
  )
}
# Para asignar las categorías tomando en cuenta los cutpoints
# y usando leave-one-identity-out
get_independ_categories <- function(DF_FEATURES, # df with all features
                                    DESVEST, # standard deviation 
                                    NOMBRES_RASGOS # features that will be categorized 
                                    ) {
  loo_results <- list()
  for (id in unique(DF_FEATURES$individual)) {
    # Separar referencia y target
    reference <- DF_FEATURES %>%
      filter(individual != id)
    
    target <- DF_FEATURES %>%
      filter(individual == id)
    # Obtener categorias de target usando la referencia
    target_cat <- target
    
    for (trait in NOMBRES_RASGOS) {
      # Calcular límites en reference
      cuts <- get_cutpoints(
        reference[[trait]],
        cutoff = DESVEST
      )
      
      # Asignar categoría a target
      target_cat[[trait]] <- case_when(
        target[[trait]] <= cuts$lower ~ "left",
        target[[trait]] >= cuts$upper ~ "right",
        TRUE ~ "middle"
      )
    }
    
    # Guardar resultados para cada identidad
    loo_results[[as.character(id)]] <- target_cat
  }
  return(loo_results)
} 


# Ajustar de manera correspondiente 
SD_CUTOFF = 1
RASGOS = TRAITS_FULL

# Asignar categorías usando leave-one-out 
features_cat <- get_independ_categories(features, SD_CUTOFF, RASGOS )
features_cat_df <- bind_rows(features_cat, .id = "omitted_id")

# CALCULO DE LR ####
## Función 
# Para calcular el LR
run_lr_analysis <- function(DF,TRAITS) {
  # Comprobar que existen los rasgos
  missing_traits <- setdiff(TRAITS, names(DF))
  
  if (length(missing_traits) > 0) {
    stop(
      "No se encontraron los siguientes rasgos: ",
      paste(missing_traits, collapse = ", ")
    )
  }
  
  # Generar todos los pares de fotografías
  pair_index <- t( combn(seq_len(nrow(DF)), 2)  )
  
  pairs <- data.frame(
    row1 = pair_index[, 1],
    row2 = pair_index[, 2],
    
    name1 = DF$name[pair_index[, 1]],
    name2 = DF$name[pair_index[, 2]],
    
    individual1 = DF$individual[pair_index[, 1]],
    individual2 = DF$individual[pair_index[, 2]]
  )
  
  # clasificar las distribuciones
  pairs$type <- ifelse(
    pairs$individual1 == pairs$individual2,
    "intra","inter")
  
  # Matriz de categorias
  trait_matrix <- as.matrix(
    DF[, TRAITS, drop = FALSE]
  )
  
  # Comparar las coincidencias por pares
  # TRUE = similitud, FALSE = disimilitud
  matches <- trait_matrix[pairs$row1, , drop = FALSE] ==
    trait_matrix[pairs$row2, , drop = FALSE]
  
  # numero acumulado de similitudes, e.g. la columna 3 indica 
  # cuántas similitudes completas hay entre los primeros 3 rasgos.
  cumulative_matches <- t(apply(matches, 1, cumsum))
  
  # P_same
  # Número de pares intraindividuales que coinciden en todos
  # los primeros k rasgos
  intra_index <- which(pairs$type == "intra")
  P_same <- sapply(
    seq_along(TRAITS),
    function(k) {
      numerator <- sum(
        # ¿El par coincide en los dos primeros k rasgos?
        cumulative_matches[intra_index, k] == k
      )
      # Número total de pares intraindividuales
      denominator <- length(intra_index)
      # Psk
      numerator / denominator
    }
  )

  # P_diff 
  # Número de pares interindividuales que coinciden en todos
  # los primeros k rasgos
  inter_index <- which(pairs$type == "inter")
  P_diff <- sapply(
    seq_along(TRAITS),
    function(k) {
      numerator <- sum(
        # ¿El par coincide en los dos primeros k rasgos?
        cumulative_matches[inter_index, k] == k
      )
      # Número total de pares intraindividuales
      denominator <- length(inter_index)
      # Pdk
      numerator / denominator
    }
  )
  # LR
  LR <- P_same / P_diff
  
  # Resultados
  results <- data.frame(
    traits = TRAITS,
    n_traits = seq_along(TRAITS),
    P_same = P_same,
    P_diff = P_diff,
    LR = LR
  )
  
  return(results)
}


# Orden predefinido: ESCENARIO CON INFORMACIÓN COMPLETA ####
results <- run_lr_analysis(
  DF = features_cat_df, TRAITS = RASGOS
  )


# Gráfica de probabilidades y LR ####
a <- ggplot(results, aes(x = n_traits, y = P_same)) +
  geom_point(size = 2, color = "#f1a226") +
  geom_point(data = results, aes(x = n_traits, y = P_diff), 
             size = 2, color ="#298c8c") +
  labs(
    x = "Number of traits",
    y = "Probability",
    title = "Complete-information scenario"
  ) + 
  theme_bw() +
  scale_x_continuous(
    breaks = 1:nrow(results) # Specify where the breaks (ticks) should be
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2) # Specify where the breaks (ticks) should be
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 32)) 

b <- ggplot(results, aes(x = n_traits, y = LR)) +
  geom_point(size = 2) + theme_bw() +
  scale_x_continuous(
    breaks = 1:nrow(results) # Specify where the breaks (ticks) should be
  )+
  scale_y_continuous(
    breaks = seq(0, 100, by = 5)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood Ratio"
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = 0.2, fill = "#0072b2")


# Orden predefinido: ESCENARIO CON INFORMACIÓN INCOMPLETA ####
# Gorra
results_simulation_cap <- run_lr_analysis(
  DF = features_cat_df,
  TRAITS = TRAITS_CAP
)
# Lentes
results_simulation_glasses <- run_lr_analysis(
  DF = features_cat_df,
  TRAITS = TRAITS_GLASSES
)
# Mascarilla
results_simulation_mask <- run_lr_analysis(
  DF = features_cat_df,
  TRAITS = TRAITS_MASK
)

# Gráfica de probabilidades y LR ####
c <- ggplot(results_simulation_cap, aes(x = n_traits, y = P_same)) +
  geom_point(size = 2, color = "#f1a226") +
  geom_point(data = results_simulation_cap, aes(x = n_traits, y = P_diff), 
             size = 2, color ="#298c8c") +
  labs(
    x = "Number of traits",
    y = "Probability",
    subtitle = "CAP"
  ) + 
  theme_bw()+
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_cap) # Specify where the breaks (ticks) should be
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2) # Specify where the breaks (ticks) should be
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 25)) 

d <- ggplot(results_simulation_cap, aes(x = n_traits, y = LR)) +
  geom_point(size = 2) + theme_bw() +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_cap) # Specify where the breaks (ticks) should be
  ) +
  scale_y_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood Ratio"
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2")

c + d

e <- ggplot(results_simulation_glasses, aes(x = n_traits, y = P_same)) +
  geom_point(size = 2, color = "#f1a226") +
  geom_point(data = results_simulation_glasses, aes(x = n_traits, y = P_diff), 
             size = 2, color ="#298c8c") +
  labs(
    x = "Number of traits",
    y = "Probability",
    title = "Reduced-information scenario",
    subtitle = "GLASSES"
  ) + 
  theme_bw()+
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_glasses) # Specify where the breaks (ticks) should be
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2) # Specify where the breaks (ticks) should be
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 25),
        plot.title = element_text(size = 32)) 

f <- ggplot(results_simulation_glasses, aes(x = n_traits, y = LR)) +
  geom_point(size = 2) + theme_bw() +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_glasses) # Specify where the breaks (ticks) should be
  )+
  scale_y_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood Ratio"
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2")

e + f

g <- ggplot(results_simulation_mask, aes(x = n_traits, y = P_same)) +
  geom_point(size = 2, color = "#f1a226") +
  geom_point(data = results_simulation_mask, aes(x = n_traits, y = P_diff), 
             size = 2, color ="#298c8c") +
  labs(
    x = "Number of traits",
    y = "Probability",
    subtitle = "MASK"
  ) + 
  theme_bw()+
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_mask) # Specify where the breaks (ticks) should be
  ) +
  scale_y_continuous(
    limits = c(0, 0.8),
    breaks = seq(0, 0.8, by = 0.2) # Specify where the breaks (ticks) should be
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 25)) 

h <- ggplot(results_simulation_mask, aes(x = n_traits, y = LR)) +
  geom_point(size = 2) + theme_bw() +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_mask) # Specify where the breaks (ticks) should be
  )+
  scale_y_continuous(
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood Ratio"
  ) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2")

g + h


FULL_LR <- ((a + b)/(e + f)/(c + d)/(g + h)) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')&
  theme(
    axis.title = element_text(size = 20),
    plot.tag = element_text(size = 25)
  )


ggsave("Figure 4.png", FULL_LR, width = 15, height = 25, units = "in", dpi = 300)


# Orden aleatorio: ESCENARIO CON INFORMACIÓN COMPLETA ####
set.seed(123)
# Permutaciones
altorder <- lapply( 1:100,
  function(i) sample(RASGOS, length(RASGOS), replace = FALSE)
)
# Calculo de LR usando las permutaciones
altorder_result <- lapply(
  seq_along(altorder),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = altorder[[i]]
    )
    
    results$order <- i
    results
  }
)

altorder_result_df <- do.call(rbind, altorder_result)
# Mediana LR
altorder_result_summary <- altorder_result_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )

# minimo y maximo por cada k
altorder_result_df %>% group_by(n_traits) %>%
  summarise(min = min(LR), 
            max = max(LR))

# Gráfica de permutaciones ####
a_alt <- ggplot(
  # LRs obtenidos en las 100 permutaciones
  altorder_result_df, aes(x = n_traits, y = LR, group = order)) +
  geom_line(alpha = 0.08, linewidth = 0.55) +
  # Mediana del LR usando las 100 permutaciones
  geom_line(data = altorder_result_summary, aes( x = n_traits, y = median_LR ),
            inherit.aes = FALSE, linewidth = 1.2 ) +
  # LRs obtenidos usando el orden predefinido
  geom_line( data = results, aes( x = n_traits, y = LR ),
             inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results)
  ) +
  scale_y_continuous(
    breaks = seq(0, 140, by = 10)
  ) +
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    title = "Complete-information scenario"
  ) +
  theme_bw() +  
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 32)
  )


# Orden aleatorio: ESCENARIO CON INFORMACIÓN INCOMPLETA ####
#### cap 
set.seed(123)
altorder_cap <- lapply(
  1:100,
  function(i) sample(TRAITS_CAP, length(TRAITS_CAP), replace = FALSE)
)
altorder_result_cap <- lapply(
  seq_along(altorder_cap),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = altorder_cap[[i]]
    )
    
    results$order <- i
    results
  }
)
altorder_result_cap_df <- do.call(rbind, altorder_result_cap)
altorder_result_summary_cap <- altorder_result_cap_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )


altorder_result_cap_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR))


#### glasses 
set.seed(123)
altorder_glasses <- lapply(
  1:100,
  function(i) sample(TRAITS_GLASSES, length(TRAITS_GLASSES), replace = FALSE)
)
altorder_result_glasses <- lapply(
  seq_along(altorder_glasses),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = altorder_glasses[[i]]
    )
    
    results$order <- i
    results
  }
)
altorder_result_glasses_df <- do.call(rbind, altorder_result_glasses)
altorder_result_summary_glasses <- altorder_result_glasses_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )

altorder_result_glasses_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR)) 


#### mask 
set.seed(123)
altorder_mask <- lapply(
  1:100,
  function(i) sample(TRAITS_MASK, length(TRAITS_MASK), replace = FALSE)
)
altorder_result_mask <- lapply(
  seq_along(altorder_mask),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = altorder_mask[[i]]
    )
    
    results$order <- i
    results
  }
)
altorder_result_mask_df <- do.call(rbind, altorder_result_mask)
altorder_result_summary_mask <- altorder_result_mask_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )

altorder_result_mask_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR))


# Gráfica de permutaciones ####
b_alt <- ggplot(altorder_result_cap_df, aes(x = n_traits, y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.55) +
  geom_line(data = altorder_result_summary_cap, aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_cap, aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_cap)
  ) +
  scale_y_continuous(
    breaks = seq(0, 35, by = 5),
    limits = c(0, 35)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    subtitle = "CAP"
  ) + 
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 25)
  )

c_alt <- ggplot(altorder_result_glasses_df,aes(x = n_traits,y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.4) +
  geom_line(data = altorder_result_summary_glasses,aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_glasses,aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_glasses)
  ) +
  scale_y_continuous(
    breaks = seq(0, 85, by = 5),
    limits = c(0, 85)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    title = "Reduced-information scenario",
    subtitle = "GLASSES"
  ) +
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 32),
    plot.subtitle = element_text(size = 25)
  )


d_alt <- ggplot(altorder_result_mask_df,aes(x = n_traits,y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.4) +
  geom_line(data = altorder_result_summary_mask,aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_mask,aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_mask)
  ) +
  scale_y_continuous(
    breaks = seq(0, 35, by = 5),
    limits = c(0, 35)
  ) +
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    subtitle = "MASK"
  ) +
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 25)
  )


design <- "AAA###
           ccbbdd"
ALTorder <- a_alt + b_alt + c_alt + d_alt + 
  plot_layout(design = design, guides = "collect") &
  theme(
    axis.title = element_text(size = 25))

ggsave("Figure 5.png", ALTorder, width = 25, height = 15, units = "in", dpi = 300)


# Selección alternativa: ESCENARIO CON INFORMACIÓN INCOMPLETA ####
#### cap
set.seed(123)
sel_cap <- lapply(
  1:100,
  function(i) sample(RASGOS, length(TRAITS_CAP), replace = FALSE)
)
sel_result_cap <- lapply(
  seq_along(sel_cap),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = sel_cap[[i]]
    )
    
    results$order <- i
    results
  }
)
sel_result_cap_df <- do.call(rbind, sel_result_cap)
sel_result_summary_cap <- sel_result_cap_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )


sel_result_cap_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR))


sel_result_cap_df[sel_result_cap_df$order == 95, ]


#### glasses ####
set.seed(123)
sel_glasses <- lapply(
  1:100,
  function(i) sample(RASGOS, length(TRAITS_GLASSES), replace = FALSE)
)
sel_result_glasses <- lapply(
  seq_along(sel_glasses),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = sel_glasses[[i]]
    )
    
    results$order <- i
    results
  }
)
sel_result_glasses_df <- do.call(rbind, sel_result_glasses)
sel_result_summary_glasses <- sel_result_glasses_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )

sel_result_glasses_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR) )


#### mask ####
set.seed(123)
sel_mask <- lapply(
  1:100,
  function(i) sample(RASGOS, length(TRAITS_MASK), replace = FALSE)
)
sel_result_mask <- lapply(
  seq_along(sel_mask),
  function(i) {
    
    results <- run_lr_analysis(
      DF = features_cat_df,
      TRAITS = sel_mask[[i]]
    )
    
    results$order <- i
    results
  }
)
sel_result_mask_df <- do.call(rbind, sel_result_mask)
sel_result_summary_mask <- sel_result_mask_df %>%
  group_by(n_traits) %>%
  summarise(
    median_LR = median(LR),
    .groups = "drop"
  )

sel_result_mask_df %>% group_by(n_traits) %>%
  summarise(min = min(LR),
            max = max(LR) )


# Gráfica de permutaciones ####
b_sel <- ggplot(sel_result_cap_df, aes(x = n_traits, y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.55) +
  geom_line(data = sel_result_summary_cap, aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_cap, aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_cap)
  ) +
  scale_y_continuous(
    breaks = seq(0, 35, by = 5),
    limits = c(0, 35)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    subtitle = "CAP"
  ) + 
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18), 
    plot.subtitle = element_text(size = 25)
  )


c_sel <- ggplot(sel_result_glasses_df,aes(x = n_traits,y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.4) +
  geom_line(data = sel_result_summary_glasses,aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_glasses,aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_glasses)
  ) +
  scale_y_continuous(
    breaks = seq(0, 85, by = 5),
    limits = c(0, 85)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    title = "Reduced-information scenario",
    subtitle = "GLASSES"
  ) +
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 32),
    plot.subtitle = element_text(size = 25)
    
  )


d_sel <- ggplot(sel_result_mask_df,aes(x = n_traits,y = LR, group = order)) +
  geom_line(alpha = 0.08,linewidth = 0.4) +
  geom_line(data = sel_result_summary_mask,aes(x = n_traits,y = median_LR),
            inherit.aes = FALSE,linewidth = 1.2) +
  geom_line(data = results_simulation_mask,aes(x = n_traits,y = LR),
            inherit.aes = FALSE,linetype = "dashed",linewidth = 1.2) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2, ymax = 10,
           alpha = .2, fill = "#0072b2") +
  scale_x_continuous(
    breaks = 1:nrow(results_simulation_mask)
  ) +
  scale_y_continuous(
    breaks = seq(0, 35, by = 5),
    limits = c(0, 35)
  )+
  labs(
    x = "Number of traits",
    y = "Likelihood ratio",
    subtitle = "MASK"
  ) +
  theme_bw() +  
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.subtitle = element_text(size = 25)
  )


design <- "ccbbdd"
ALTselect <- b_sel + c_sel + d_sel + 
  plot_layout(design = design, guides = "collect") &
  theme(axis.title = element_text(size = 25))

ggsave("Figure 6.png", ALTselect, width = 20, height = 8, units = "in", dpi = 300)

