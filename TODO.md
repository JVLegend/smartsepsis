# TODO - SmartLab BacEnd

Tarefas pendentes e em andamento do projeto.
Ultima atualizacao: 2026-04-16

---

## Concluido

- [x] Pipeline base: fetch, guides, primers, BLAST, panel (5 etapas)
- [x] 12 familias genicas de referencia processadas com BLAST 100%
- [x] 30 variantes clinicas processadas
- [x] Conservation analysis para cobertura de variantes
- [x] Sistema de tracking CSV
- [x] Batch mode (run_batch.py)
- [x] Site vitrine (public/index.html)
- [x] Documentacao: README.md, TODO.md, MELHORIAS.md

---

## Em andamento

### Fase 1 - Integracao CARD Database [CONCLUIDA]
- [x] Criar modulo `card_integration.py` para buscar dados do CARD (Comprehensive Antibiotic Resistance Database)
- [x] Mapear variantes CARD -> alvos existentes do pipeline (12 familias)
- [x] Auto-descoberta de novas variantes AMR relevantes para o Brasil
- [x] Enriquecer dados com metadados CARD (mecanismo de resistencia, drug class, ontologia ARO)
- [x] Gerar CSV enriquecido (`targets_brazil_card.csv`) e relatorio de descobertas
- [x] Atualizar site com secao CARD

### Fase 2 - Explicacoes em Linguagem Natural (Claude API) [CONCLUIDA]
- [x] Criar modulo `clinical_interpreter.py` com integracao Claude API
- [x] Gerar interpretacao clinica para cada alvo/guide/primer (3 audiencias: medico, gestor, pesquisador)
- [x] Contexto clinico embutido para 12 familias AMR (mortalidade, isolamento, tratamento, prevalencia BR)
- [x] Modo online (Claude API) + modo offline (pre-computado)
- [x] Relatorio clinico em .md + JSON para frontend
- [x] Atualizar site com secao "Interpretacao Clinica por IA"

### Fase 3 - Score Funcional via Evo 2 [CONCLUIDA]
- [x] Criar modulo `evo2_scoring.py` com modo GPU (Evo 2 real) e lightweight (features de sequencia)
- [x] 5 dimensoes de analise: k-mer profile, codon adaptation, GC shift, dinucleotideo, complexidade
- [x] Calcular distancia funcional composta entre variantes e referencia
- [x] Classificacao de impacto: conserved, likely_conserved, uncertain, likely_disrupted, disrupted
- [x] Priorizacao dinamica baseada em impacto funcional
- [x] Outputs: TSV, JSON, relatorio textual
- [x] Atualizar site com secao "Score Funcional de Variantes"

### Fase 4 - Dashboard Interativo [CONCLUIDA]
- [x] Criar `public/dashboard.html` - pagina interativa completa
- [x] Tabela sortavel com 8 colunas (gene, classe, prioridade, organismo, variantes, BLAST, impacto, acao)
- [x] Filtros por prioridade (P1/P2/P3), classe de resistencia e busca por texto
- [x] Painel de detalhes expandivel com 3 colunas: interpretacao clinica, dados pipeline, dados CARD
- [x] Mapa visual do gene mostrando posicao do guide e primers
- [x] Barra de impacto funcional por alvo
- [x] Exportacao CSV
- [x] Botao "Abrir Dashboard" no hero do site principal
- [x] Navegacao entre dashboard e site principal

### Fase 5 - Covariance Probes para Design de Guias [CONCLUIDA]
- [x] Criar modulo `covariance_probes.py` com 18 features biofisicas
- [x] Features novas: seed region (GC, purinas, homopolymer), termodinamica (dG nearest-neighbor, energia seed, gradiente estabilidade), acessibilidade estrutural (palindromes, repeticoes invertidas), preferencias posicionais (dados experimentais Cas12a)
- [x] Implementar covariance matrix sobre todos os candidatos
- [x] Score composto: w^T*f + alpha * f^T*C*w (linear + covariancia)
- [x] Benchmark automatico contra scoring rule-based original
- [x] Outputs: TSV, JSON, relatorio de benchmark detalhado
- [x] Atualizar site com secao "Covariance Probes" + pipeline atualizado para 7 etapas

### Fase 6 - 15 Analises Complementares [CONCLUIDA - Abr/2026]
- [x] Guide redesign termodinamica (seed+PAM) - 42/42 guides com GC 40-60%
- [x] Embedding comparison (Evo 2 distances) - heatmap
- [x] Dashboard interativo de conservacao - HTML com charts
- [x] Cross-reactivity analysis - 0 off-target hits >=18nt
- [x] Publication stats (LaTeX tables + CI)
- [x] Codon usage bias (RSCU) - rscu_analysis.csv (72KB)
- [x] Hamming distance matrix - hamming_distance_matrix.csv (13KB)
- [x] Coverage heatmap gene-family x match-type
- [x] GC/Tm biofisica - Tm 54-67C (RPA-compativel)
- [x] CARD enrichment - card_enrichment.html (111KB) - descobriu KPC-33, NDM-9, OXA-244
- [x] Comprehensive multi-secao report
- [x] Multiplex panel optimization - 12 guides cobrindo 42 alvos
- [x] Network adjacency gene-gene
- [x] Diversity statistics (Shannon, nucleotide diversity)
- [x] Summary dashboard executivo (KPIs)
- [x] Atualizar site com secao "Melhorias Experimentais e Resultados"
- [x] Criar AUTORESEARCH_DGX_RESULTADOS.md

---

## Em andamento / Bloqueado

### Copia de resultados DGX -> Mac
- [ ] Baixar 15 outputs HTML/CSV quando SSH voltar
  ```bash
  scp -P 2222 oftalmousp@dgx.retina.ia.br:~/jv-teste/harvard_bacend/results/*.{csv,html} \
      /Users/jv/Documents/GitHub/bac_end_harvard/results_dgx/
  ```

---

## Proximas iteracoes (Fase 7+)

- [ ] Incluir KPC-33, NDM-9, OXA-244 no painel
- [ ] Validacao experimental de cross-reactivity em cepas HC-FMUSP
- [ ] Phylogenetic tree por familia (ClustalW + iqtree)
- [ ] ESM-2 embeddings (1280d) para cada proteina AMR
- [ ] Tajima's D test de selecao por gene
- [ ] dN/dS ratio analysis
- [ ] HGT (horizontal gene transfer) detection via codon incongruence
- [ ] ILP multiplex optimization formal
- [ ] RNAfold secondary structure analysis (evitar hairpins)
- [ ] Regional-specific variants analysis (Brasil vs mundo)

---

## Backlog (futuro)

- [ ] Integracao com dados epidemiologicos em tempo real (ANVISA/BR-GLASS API)
- [ ] Suporte a novos patogenos alem de AMR (ex: virus respiratorios)
- [ ] App mobile para leitura de fluorescencia via camera
- [ ] Validacao clinica com amostras reais (parceria HC-FMUSP)
- [ ] Submissao regulatoria ANVISA (IVD, RDC 830/2023)
- [ ] Multi-idioma (ingles/espanhol) para adocao internacional
