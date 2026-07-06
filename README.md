# DIABLO – Integração Multiômica Supervisionada

Este repositório contém um exemplo de pipeline para integração de dados multiômicos utilizando o método **DIABLO (Data Integration Analysis for Biomarker discovery using Latent Components)**, implementado no pacote **mixOmics**.

O objetivo do script é demonstrar o fluxo básico de uma análise supervisionada envolvendo transcriptômica, proteômica e metilação de DNA para identificação de assinaturas moleculares associadas a um fenótipo de interesse.

## Fluxo da análise

O pipeline contempla as seguintes etapas:

1. Construção dos blocos de dados (`X`) e variável resposta (`Y`);
2. Definição da matriz de desenho (*design matrix*);
3. Definição do espaço de busca para seleção de variáveis (`keepX`);
4. Otimização do modelo por validação cruzada;
5. Ajuste do modelo final (`block.splsda`);
6. Visualização dos resultados:
   - projeção das amostras;
   - correlação entre blocos;
   - arrow plot;
   - correlation circle plot;
   - circos plot;
   - loadings;
7. Avaliação do desempenho por curvas ROC e AUC;
8. Avaliação da taxa de erro por validação cruzada.

## Dependências

O pipeline utiliza principalmente os seguintes pacotes:

- mixOmics
- pROC
- ggplot2
- dplyr

## Observações

Este script foi desenvolvido como material de apoio para uma disciplina de análise multiômica e tem como objetivo ilustrar a aplicação do método DIABLO. Etapas de pré-processamento dos dados (filtragem, normalização, seleção de amostras e controle de qualidade) não estão incluídas neste repositório, uma vez que dependem das características de cada conjunto de dados.

---
**Autor:** Luiz Garcia
