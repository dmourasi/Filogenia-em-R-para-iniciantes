# Filogenia para Iniciantes em R

Este tutorial foi criado por Olivier De Clerk (Ugent). Atualizado e expandido por Daniel Moura.

---

## Sumário

1. [Introdução à Filogenética](#1-introdução-à-filogenética)
2. [Instalando os pacotes necessários](#2-instalando-os-pacotes-necessários)
3. [Árvores: estruturas básicas com `ape`](#3-árvores-estruturas-básicas-com-ape)
4. [Formatos de arquivo: Newick e Nexus](#4-formatos-de-arquivo-newick-e-nexus)
5. [Escrevendo e lendo árvores filogenéticas](#5-escrevendo-e-lendo-árvores-filogenéticas)
6. [Alinhamento de sequências](#6-alinhamento-de-sequências)
7. [Métodos baseados em distância: UPGMA e Neighbor-Joining](#7-métodos-baseados-em-distância-upgma-e-neighbor-joining)
8. [Máxima Parcimônia](#8-máxima-parcimônia)
9. [Máxima Verossimilhança com `phangorn`](#9-máxima-verossimilhança-com-phangorn)
10. [Suporte de ramos: Bootstrap](#10-suporte-de-ramos-bootstrap)
11. [Inferência Bayesiana](#11-inferência-bayesiana)
12. [Visualização avançada de árvores](#12-visualização-avançada-de-árvores)
13. [Reconstrução de estados ancestrais](#13-reconstrução-de-estados-ancestrais)
14. [Métodos comparativos filogenéticos](#14-métodos-comparativos-filogenéticos)
15. [Referências e leituras adicionais](#15-referências-e-leituras-adicionais)

---

## 1. Introdução à Filogenética

A **filogenética** é o estudo das relações evolutivas entre organismos. Uma **árvore filogenética** (ou filogenia) é uma representação gráfica dessas relações, onde:

- **Folhas (tips/OTUs)**: representam os táxons atuais (espécies, populações, sequências)
- **Nós internos**: representam ancestrais comuns hipotéticos
- **Ramos (branches/edges)**: conectam nós e podem ter comprimentos que representam tempo evolutivo ou número de mudanças
- **Raiz**: o ancestral comum mais antigo de todos os táxons analisados

### Árvores enraizadas vs. não-enraizadas

Uma árvore **enraizada** possui um ponto de referência (a raiz), indicando a direção do tempo evolutivo. Uma árvore **não-enraizada** mostra apenas as relações entre os táxons sem indicar qual grupo é mais basal. Para enraizar uma árvore, geralmente usamos um **grupo externo (outgroup)** — um táxon sabidamente mais distante de todos os outros no estudo.

### Principais métodos de inferência filogenética

| Método | Critério | Pacote R |
|--------|----------|----------|
| UPGMA / Neighbor-Joining | Distância | `ape` |
| Máxima Parcimônia | Parcimônia | `phangorn` |
| Máxima Verossimilhança | Verossimilhança | `phangorn`, `RAxML` |
| Inferência Bayesiana | Probabilidade posterior | `MrBayes`, `BEAST` |

---

## 2. Instalando os pacotes necessários

R oferece uma ampla gama de pacotes para filogenética. Uma visão geral quase completa pode ser encontrada em: https://cran.r-project.org/web/views/Phylogenetics.html

Neste tutorial usaremos principalmente:

```r
# Instale todos os pacotes de uma vez (execute apenas uma vez por máquina)
install.packages(c("ape", "phytools", "phangorn", "ggtree", "ggplot2", "msa", "seqinr"))

# Para instalar o pacote ggtree (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("msa")
```

Carregar as bibliotecas no início de cada sessão:

```r
library(ape)
library(phytools)
library(phangorn)
library(ggplot2)
```

---

## 3. Árvores: estruturas básicas com `ape`

O pacote `ape` (*Analyses of Phylogenetics and Evolution*) é o núcleo da filogenética em R. A maioria dos outros pacotes depende das estruturas de dados definidas no `ape`.

> Referência principal: Paradis, E. 2011. *Analysis of Phylogenetics and Evolution with R*. Springer Science & Business Media.

### 3.1 Simulando e visualizando árvores

```r
library(ape)

# Simula uma árvore aleatória com 20 táxons
tree <- rtree(n = 20)
plot(tree, edge.width = 2)
```

![image](https://user-images.githubusercontent.com/23074735/89691809-4a5b6900-d8e0-11ea-8218-423dafa12dc8.png)

A função `rtree()` gera uma árvore aleatória com topologia e comprimentos de ramos aleatórios. O argumento `n` define o número de táxons (folhas).

### 3.2 A estrutura interna de um objeto `phylo`

Um objeto da classe `phylo` (o formato padrão do `ape`) é uma **lista** com pelo menos 4 componentes. Vamos explorar isso com um exemplo simples:

```r
# Lê uma árvore no formato Newick a partir de texto
tree <- read.tree(text = "(((A,B),(C,D)),E);")
plot(tree, edge.width = 2, type = "cladogram", font = 4)
```

Inspecionando a estrutura interna:

```r
# Matriz de arestas: cada linha é [nó_ancestral, nó_descendente]
tree$edge

# Rótulos das folhas (na mesma ordem das pontas 1 a n)
tree$tip.label

# Comprimentos dos ramos (se existirem)
tree$edge.length

# Número de táxons
tree$Nnode
```

Visualizando a numeração de nós e pontas:

```r
plot(tree, type = "phylogram", edge.width = 2)
tiplabels()   # numera as pontas
nodelabels()  # numera os nós internos
```

**Convenção de numeração do `ape`:**
- Pontas (folhas): numeradas de `1` a `n`
- Nós internos: numerados de `n+1` a `n+m`, onde `m = n-1` para árvores bifurcadas completas
- A raiz é sempre o nó `n+1`

### 3.3 Tipos de visualização

```r
par(mfrow = c(2, 3))  # divide a janela em 6 painéis

tree <- rtree(15)
plot(tree, type = "phylogram",  main = "phylogram")
plot(tree, type = "cladogram",  main = "cladogram")
plot(tree, type = "fan",        main = "fan")
plot(tree, type = "unrooted",   main = "unrooted")
plot(tree, type = "radial",     main = "radial")

par(mfrow = c(1, 1))  # restaura o painel único
```

---

## 4. Formatos de arquivo: Newick e Nexus

Antes de trabalhar com dados reais, é fundamental entender os formatos de arquivo mais comuns.

### 4.1 Formato Newick

O formato **Newick** (também chamado de New Hampshire) representa uma árvore como uma sequência de parênteses aninhados. É o formato mais simples e amplamente suportado.

```
# Árvore sem comprimentos de ramo
(((A,B),(C,D)),E);

# Árvore com comprimentos de ramo
(((A:0.1,B:0.2):0.05,(C:0.3,D:0.15):0.08):0.12,E:0.4);

# Árvore com suporte de ramos e comprimentos
(((A:0.1,B:0.2)90:0.05,(C:0.3,D:0.15)85:0.08)95:0.12,E:0.4);
```

Lendo e escrevendo no formato Newick:

```r
# Lendo de um arquivo
tree <- read.tree("minha_arvore.nwk")

# Lendo de uma string de texto
tree <- read.tree(text = "(((A:0.1,B:0.2),C:0.3),D:0.5);")

# Escrevendo em arquivo
write.tree(tree, "minha_arvore.nwk")

# Imprimindo a string Newick no console
write.tree(tree)
```

### 4.2 Formato Nexus

O formato **Nexus** é mais rico e permite armazenar não apenas árvores, mas também alinhamentos, matrizes de caracteres e blocos de comandos. É o formato padrão para programas como MrBayes e PAUP*.

```
#NEXUS

BEGIN TAXA;
    NTAX = 4;
    TAXLABELS A B C D;
END;

BEGIN TREES;
    TREE 1 = (((A:0.1,B:0.2):0.05,C:0.3):0.08,D:0.5);
END;
```

Lendo e escrevendo no formato Nexus:

```r
library(phytools)

# Lendo
tree <- read.nexus("minha_arvore.nex")

# Escrevendo
writeNexus(tree, "minha_arvore.nex")

# Visualizando o conteúdo do arquivo
cat(readLines("minha_arvore.nex"), sep = "\n")
```

---

## 5. Escrevendo e lendo árvores filogenéticas

### 5.1 Operações básicas com `phytools`

```r
library(phytools)

# Simulando uma árvore com modelo birth-death
# b = taxa de especiação, d = taxa de extinção
tree <- pbtree(b = 1, d = 0.2, n = 40)
plotTree(tree)
nodelabels()
```

### 5.2 Extraindo clados e podando pontas

```r
# Extraindo o clado descendente de um nó específico
# (verifique a numeração dos nós com nodelabels() antes)
clado <- extract.clade(tree, node = 60)
plotTree(clado)

# Salvando o clado como arquivo separado
writeNexus(clado, "meu_clado.nex")
```

```r
# Removendo pontas aleatórias da árvore
dtips <- sample(tree$tip.label, 10)
dt <- drop.tip(tree, dtips)
plotTree(dt)
```

```r
# Removendo apenas as pontas extintas (linhagens que não chegam ao presente)
et <- fancyTree(tree, type = "droptip", tip = getExtinct(tree), cex = 0.7)
```

### 5.3 Manipulações úteis de árvores

```r
# Enraizar uma árvore não-enraizada com um outgroup
tree_unrooted <- read.tree(text = "((A,B),(C,(D,E)));")
tree_rooted <- root(tree_unrooted, outgroup = "A", resolve.root = TRUE)

# Descartar a raiz (tornar não-enraizada)
tree_unr <- unroot(tree_rooted)

# Ladderizar (ordenar os ramos para melhor visualização)
tree_lad <- ladderize(tree_rooted)
plot(tree_lad)

# Obter todos os ancestrais de uma ponta
getMRCA(tree_rooted, tip = c("C", "D"))  # Ancestral comum mais recente

# Obter a distância patrística entre pontas
cophenetic(tree_rooted)
```

---

## 6. Alinhamento de sequências

Na prática, a análise filogenética começa com **sequências de DNA, RNA ou proteínas**. Antes de inferir a filogenia, as sequências precisam ser **alinhadas** para identificar posições homólogas.

### 6.1 Lendo sequências no formato FASTA

```r
library(ape)

# Lendo sequências de DNA no formato FASTA
seqs <- read.dna("minhas_sequencias.fasta", format = "fasta")

# Verificando o objeto
class(seqs)    # "DNAbin"
dim(seqs)      # número de sequências x comprimento
rownames(seqs) # nomes dos táxons
```

### 6.2 Realizando o alinhamento múltiplo

```r
library(msa)

# Lendo sequências não-alinhadas
seqs_raw <- readDNAStringSet("sequencias_nao_alinhadas.fasta")

# Alinhamento com ClustalW
aln_clustal <- msa(seqs_raw, method = "ClustalW")

# Alinhamento com MUSCLE (geralmente mais preciso)
aln_muscle <- msa(seqs_raw, method = "Muscle")

# Convertendo para formato DNAbin do ape
aln_ape <- as.DNAbin(aln_muscle)
```

### 6.3 Convertendo para o formato `phyDat`

O pacote `phangorn` usa o formato `phyDat` para os dados de sequência:

```r
library(phangorn)

# Convertendo DNAbin para phyDat
dados <- phyDat(aln_ape, type = "DNA")

# Para proteínas
# dados <- phyDat(aln_ape, type = "AA")
```

---

## 7. Métodos baseados em distância: UPGMA e Neighbor-Joining

Os métodos de distância transformam as sequências em uma **matriz de distâncias** e depois constroem a árvore a partir dessa matriz. São computacionalmente rápidos e úteis para uma análise exploratória inicial.

### 7.1 Calculando a matriz de distâncias

```r
library(ape)
library(phangorn)

# Usando dados do phangorn (exemplo embutido — mamíferos placentários)
data("Laurasiatherian", package = "phangorn")
dados <- Laurasiatherian
aln   <- as.DNAbin(dados)

# Calculando distâncias com o modelo de Jukes-Cantor (JC69)
dist_jc <- dist.dna(aln, model = "JC69")

# Outros modelos disponíveis: "K80", "F84", "HKY85", "T92", "TN93", "GG95", "logdet", "raw"
dist_k80 <- dist.dna(aln, model = "K80")
```

### 7.2 UPGMA

O **UPGMA** (*Unweighted Pair Group Method with Arithmetic mean*) assume que todos os táxons evoluem na mesma taxa (relógio molecular estrito). É simples, mas essa suposição raramente é verdadeira em dados reais.

```r
# Construindo árvore UPGMA
tree_upgma <- upgma(dist_jc)
plot(tree_upgma, main = "UPGMA")
```

### 7.3 Neighbor-Joining

O **Neighbor-Joining (NJ)** é um método muito mais robusto que não assume taxa igual de evolução. É amplamente usado como método rápido de análise exploratória.

```r
# Construindo árvore Neighbor-Joining
tree_nj <- nj(dist_jc)
plot(tree_nj, main = "Neighbor-Joining")

# Comparando UPGMA e NJ
par(mfrow = c(1, 2))
plot(tree_upgma, main = "UPGMA")
plot(tree_nj,    main = "Neighbor-Joining")
par(mfrow = c(1, 1))
```

---

## 8. Máxima Parcimônia

O critério de **Máxima Parcimônia (MP)** busca a árvore que requer o menor número de mudanças evolutivas para explicar os dados observados — o princípio da navalha de Occam aplicado à evolução.

```r
library(phangorn)

# Carregando dados de exemplo
data("Laurasiatherian", package = "phangorn")
dados <- Laurasiatherian

# Gerando uma árvore inicial com NJ para usar como ponto de partida
dist_mat   <- dist.hamming(dados)
tree_start <- NJ(dist_mat)

# Calculando o score de parcimônia de uma árvore
parsimony(tree_start, dados)

# Busca pela árvore de máxima parcimônia
# método "optim.parsimony" faz busca local (rearranjos SPR/NNI)
tree_mp <- optim.parsimony(tree_start, dados)
parsimony(tree_mp, dados)

# Busca mais exaustiva com pratchet (perturbação + otimização)
tree_mp2 <- pratchet(dados, start = tree_start, maxit = 100, minit = 20)

# pratchet pode retornar múltiplas árvores; obter consenso se necessário
if (inherits(tree_mp2, "multiPhylo")) tree_mp2 <- consensus(tree_mp2, p = 0.5)
parsimony(tree_mp2, dados)

plot(tree_mp2, main = "Máxima Parcimônia")
```

---

## 9. Máxima Verossimilhança com `phangorn`

A **Máxima Verossimilhança (ML)** é o padrão ouro entre os métodos não-Bayesianos. Ela busca a árvore (topologia + comprimentos de ramo) que maximiza a probabilidade de observar os dados sob um modelo de substituição nucleotídica explícito.

### 9.1 Seleção do modelo de substituição

Antes de inferir a árvore, é importante identificar o modelo que melhor descreve os dados usando critérios de informação (AIC, BIC):

```r
library(phangorn)

data("Laurasiatherian", package = "phangorn")
dados    <- Laurasiatherian
tree_ini <- NJ(dist.hamming(dados))

# Testando todos os modelos e selecionando pelo AIC
mt <- modelTest(dados, tree = tree_ini, model = "all")
print(mt)

# O melhor modelo (menor AIC) pode ser extraído:
melhor_modelo <- mt$Model[which.min(mt$AIC)]
cat("Melhor modelo:", melhor_modelo, "\n")
```

Modelos comuns de substituição de DNA:

| Modelo | Parâmetros | Descrição |
|--------|-----------|-----------|
| JC69   | 1 | Jukes-Cantor: todas as taxas iguais |
| K80    | 2 | Kimura 2-parâmetros: ts ≠ tv |
| HKY85  | 4 | Freqüências de base desiguais + ts ≠ tv |
| GTR    | 9 | General Time Reversible: modelo mais geral |
| GTR+G  | 10 | GTR + variação de taxa entre sítios (distribuição Gamma) |
| GTR+G+I| 11 | GTR + Gamma + proporção de sítios invariáveis |

### 9.2 Inferindo a árvore por ML

```r
# Criando o objeto pml com a árvore inicial e o modelo
fit_ini <- pml(tree_ini, dados)

# Otimizando com o modelo GTR+G+I
fit_ml <- optim.pml(fit_ini,
                    model    = "GTR",
                    optGamma = TRUE,    # otimiza variação de taxa (parâmetro alpha)
                    optInv   = TRUE,    # otimiza proporção de sítios invariáveis
                    rearrangement = "stochastic",  # busca mais agressiva de topologia
                    control  = pml.control(trace = 0))

# Resumo do ajuste
summary(fit_ml)

# Visualizando a árvore ML
plot(fit_ml$tree, main = "Máxima Verossimilhança (GTR+G+I)")
add.scale.bar()
```

---

## 10. Suporte de Ramos: Bootstrap

O **bootstrap não-paramétrico** é a medida de suporte mais usada em filogenética. Ele reamostras as colunas do alinhamento com reposição e reconstrói a árvore muitas vezes. O **valor de bootstrap** de um ramo indica em quantas porcentagem das réplicas aquele ramo apareceu.

```r
library(phangorn)

# Bootstrap com 100 réplicas (use 1000 para publicação)
bs <- bootstrap.pml(fit_ml,
                    bs      = 100,
                    optNni  = TRUE,
                    control = pml.control(trace = 0))

# Plotando a árvore consenso com valores de bootstrap
tree_bs <- plotBS(fit_ml$tree, bs, type = "phylogram", bs.col = "blue")

# Também é possível plotar a árvore ML com os valores de bootstrap
plotBS(midpoint(fit_ml$tree), bs, type = "fan")
```

**Interpretação:** Valores de bootstrap ≥ 70% são geralmente considerados como suporte moderado a forte; valores ≥ 90% indicam suporte robusto.

---

## 11. Inferência Bayesiana

A **Inferência Bayesiana (IB)** estima a **distribuição de probabilidade posterior** das topologias e parâmetros do modelo usando o Teorema de Bayes e amostragem via **MCMC** (Markov Chain Monte Carlo). O resultado é um conjunto de árvores cuja frequência representa a probabilidade posterior de cada topologia.

### 11.1 Visão geral do método

$$P(\text{árvore} | \text{dados}) \propto P(\text{dados} | \text{árvore}) \times P(\text{árvore})$$

- **Verossimilhança** `P(dados | árvore)`: calculada com um modelo de substituição
- **Prior** `P(árvore)`: conhecimento ou suposições anteriores sobre a árvore
- **Posterior** `P(árvore | dados)`: o que queremos estimar

### 11.2 MrBayes (externo ao R)

**MrBayes** é o programa de IB mais popular para filogenética molecular. É executado externamente ao R, mas o arquivo de entrada (`.nex`) pode ser preparado pelo R:

```r
library(ape)

# Escrevendo o arquivo Nexus com o bloco MrBayes embutido
writeNexus(aln_ape, "analise_mrbayes.nex")

# Adicionar manualmente ao arquivo o bloco:
# BEGIN MRBAYES;
#   set autoclose=yes;
#   lset nst=6 rates=invgamma;   # modelo GTR+G+I
#   mcmc ngen=1000000 samplefreq=500 printfreq=500;
#   sump burnin=250;
#   sumt burnin=250;
# END;
```

### 11.3 Lendo os resultados do MrBayes no R

```r
library(ape)

# Lendo a árvore consenso gerada pelo MrBayes
tree_bayes <- read.nexus("analise_mrbayes.nex.con.tre")
plot(tree_bayes, main = "Inferência Bayesiana (MrBayes)")

# Os valores nos nós são probabilidades posteriores (0 a 1)
nodelabels(tree_bayes$node.label, cex = 0.7)
```

### 11.4 BEAST2 (para análises tempo-calibradas)

**BEAST2** é o padrão para análises **filogeográficas** e **filogenias calibradas por tempo** (árvores ultramétricas com datas absolutas). Ele também roda externamente e seus resultados podem ser importados para R com o pacote `treeio` ou visualizados com `ggtree`.

---

## 12. Visualização Avançada de Árvores

### 12.1 Personalização com `ape`

```r
library(ape)

tree <- rtree(20)
tree <- ladderize(tree)

# Cores e espessuras por clado
is_clade <- c(rep("red", 10), rep("blue", 10))  # exemplo simples

plot(tree,
     tip.color = is_clade,
     edge.width = 2,
     edge.color = "darkgray",
     label.offset = 0.02,
     cex = 0.8,
     main = "Árvore personalizada")

# Adicionando escala de tempo
add.scale.bar(cex = 0.8)

# Destacando um clado com um retângulo
rect(0, 0.5, 1, 10.5, col = rgb(0.5, 0.8, 0.5, 0.3), border = NA)
```

### 12.2 Visualização com `ggtree` (baseado em `ggplot2`)

O pacote `ggtree` oferece uma interface baseada em `ggplot2` para criar visualizações de alta qualidade:

```r
library(ggtree)
library(ggplot2)

tree <- rtree(30)

# Árvore básica
ggtree(tree) + geom_tiplab(size = 3)

# Árvore circular
ggtree(tree, layout = "circular") + geom_tiplab(size = 2)

# Adicionando metadados
# Criando uma tabela de dados associada às pontas
metadata <- data.frame(
  label  = tree$tip.label,
  grupo  = sample(c("A", "B", "C"), 30, replace = TRUE),
  valor  = runif(30)
)

# Plotando com cores por grupo
ggtree(tree) %<+% metadata +
  geom_tiplab(aes(color = grupo), size = 3) +
  geom_tippoint(aes(color = grupo, size = valor)) +
  scale_color_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) +
  theme_tree2() +
  labs(title = "Árvore filogenética com metadados")
```

### 12.3 Anotando clados com `phytools`

```r
library(phytools)

tree <- pbtree(n = 40)
tree <- ladderize(tree)

plotTree(tree, ftype = "i", fsize = 0.6)

# Destacando um clado específico
node <- findMRCA(tree, c("t1", "t5"))  # ajuste os nomes
cladelabels(tree, node = node, label = "Clado A", offset = 2)
```

---

## 13. Reconstrução de Estados Ancestrais

A **reconstrução de estados ancestrais** (ancestral state reconstruction) estima os valores de características fenotípicas ou moleculares em nós ancestrais da árvore.

### 13.1 Caracteres discretos

```r
library(ape)
library(phytools)

# Simulando uma árvore e um caráter discreto
tree   <- pbtree(n = 50, scale = 1)
x      <- fastBM(tree, sig2 = 1)
x_disc <- setNames(ifelse(x > 0, "A", "B"), names(x))

# Reconstrução por ML com modelo Mk (Equal Rates) via ace()
anc_ace <- ace(x_disc, tree, type = "discrete", model = "ER")

# Visualizando: tortas nos nós representam probabilidade posterior de cada estado
plot(tree, type = "phylogram", cex = 0.5, label.offset = 0.05)
nodelabels(pie    = anc_ace$lik.anc,
           piecol = c("firebrick", "steelblue"),
           cex    = 0.5)
legend("topleft", legend = c("Estado A", "Estado B"),
       fill = c("firebrick", "steelblue"), bty = "n")
```

### 13.2 Caracteres contínuos

```r
library(phytools)

# Simulando um caráter contínuo com Movimento Browniano
tree <- pbtree(n = 40, scale = 1)
x    <- fastBM(tree, sig2 = 0.5)

# Reconstrução por máxima verossimilhança (Movimento Browniano)
fit <- fastAnc(tree, x, vars = TRUE, CI = TRUE)

# Visualizando com escala de cor (phenogram)
phenogram(tree, x, spread.labels = TRUE, spread.cost = c(1, 0))

# Projetando o caráter na árvore com gradiente de cor
contMap(tree, x, plot = TRUE)
```

---

## 14. Métodos Comparativos Filogenéticos

Os **métodos comparativos filogenéticos** (PCMs — *Phylogenetic Comparative Methods*) controlam a não-independência estatística entre espécies causada pelo compartilhamento de histórico evolutivo. Dois táxons próximos tendem a ser mais similares simplesmente por serem aparentados, e não necessariamente por causas ecológicas.

### 14.1 Sinal filogenético

O **sinal filogenético** mede o quanto uma característica segue o padrão esperado pela filogenia. O **K de Blomberg** é a métrica mais usada: K = 1 indica evolução por Movimento Browniano; K > 1 indica conservadorismo filogenético; K < 1 indica evolução convergente.

```r
library(phytools)

tree <- pbtree(n = 50, scale = 1)
x    <- fastBM(tree)  # simula um caráter com Movimento Browniano

# Testando sinal filogenético com permutação
k_test <- phylosig(tree, x, method = "K", test = TRUE, nsim = 1000)
print(k_test)  # K e valor de p

# Lambda de Pagel: outro teste de sinal filogenético
lambda_test <- phylosig(tree, x, method = "lambda", test = TRUE)
print(lambda_test)
```

### 14.2 PGLS: Regressão com controle filogenético

O **PGLS** (*Phylogenetic Generalized Least Squares*) é a versão filogenética da regressão linear:

```r
library(ape)
library(nlme)

# Simulando dados
tree <- rtree(50)
x    <- fastBM(tree)
y    <- x * 0.8 + rnorm(50, 0, 0.5)

# Criando a estrutura de correlação filogenética (Movimento Browniano)
cor_bm <- corBrownian(phy = tree)

# PGLS com gls do pacote nlme
dados_df <- data.frame(y = y, x = x, row.names = tree$tip.label)
modelo_pgls <- gls(y ~ x, correlation = cor_bm, data = dados_df)
summary(modelo_pgls)
```

### 14.3 Modelos de evolução: BM, OU, EB

```r
library(phytools)

tree <- pbtree(n = 50, scale = 1)
x    <- fastBM(tree)

# Ajustando modelos de evolução com fitContinuous (pacote geiger)
library(geiger)

fit_bm <- fitContinuous(tree, x, model = "BM")   # Movimento Browniano
fit_ou <- fitContinuous(tree, x, model = "OU")   # Ornstein-Uhlenbeck (seleção estabilizadora)
fit_eb <- fitContinuous(tree, x, model = "EB")   # Early Burst (radiação adaptativa)

# Comparando modelos pelo AIC
# fitContinuous retorna objetos da classe "gfit"; o AIC está em $opt$aic
fit_bm$opt$aic
fit_ou$opt$aic
fit_eb$opt$aic
```

---

## 15. Referências e Leituras Adicionais

### Livros

- **Paradis, E. (2011).** *Analysis of Phylogenetics and Evolution with R* (2ª ed.). Springer. — Referência principal para o pacote `ape`.
- **Felsenstein, J. (2004).** *Inferring Phylogenies*. Sinauer Associates. — Fundamentos teóricos clássicos.
- **Yang, Z. (2014).** *Molecular Evolution: A Statistical Approach*. Oxford University Press. — Modelos de substituição e ML.
- **Revell, L.J. & Harmon, L.G. (2022).** *Phylogenetic Comparative Methods in R*. Princeton University Press.

### Tutoriais online

- [CRAN Task View: Phylogenetics](https://cran.r-project.org/web/views/Phylogenetics.html)
- [Phytools Blog (Liam Revell)](http://blog.phytools.org/)
- [Documentação do ggtree](https://yulab-smu.top/treedata-book/)
- [Tutorial do phangorn](https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html)

### Pacotes R mencionados

| Pacote | Função principal |
|--------|-----------------|
| `ape` | Estruturas de dados e funções fundamentais |
| `phytools` | Visualização e métodos comparativos |
| `phangorn` | ML, parcimônia, modelos de substituição |
| `ggtree` | Visualização avançada baseada em ggplot2 |
| `geiger` | Modelos de diversificação e evolução de caracteres |
| `nlme` | PGLS e modelos mistos |
| `msa` | Alinhamento múltiplo de sequências |

---

*Este material é destinado a estudantes de graduação e pós-graduação em biologia, bioinformática e áreas correlatas. Sugestões e correções são bem-vindas via Issues no repositório.*
