# Harvard BacEnd — Autoresearch DGX (Gemma4)

**Data**: 15-16 Abril 2026
**Infraestrutura**: DGX Spark (NVIDIA GB10, 119GB RAM) + Gemma4 via Ollama
**Status**: **15/15 tasks COMPLETAS**

#Harvard #CRISPR #BacEnd #AMR #PhD #Tecnologia #HCFMUSP #Academia

> [!note] Pipeline CRISPR-Cas12a para diagnostico de resistencia antimicrobiana (AMR). Autoresearch Gemma4 rodou 15 analises complementares em sequencia.

---

## Contexto

Projeto **BacEnd (Harvard Hackathon 2026)**: sistema de deteccao rapida de genes de resistencia antimicrobiana (AMR) usando CRISPR-Cas12a como biosensor. O pipeline DGX executou 15 tarefas de analise bioinformatica complementar via Gemma4.

### Dataset base
- **12 familias de genes AMR**: mecA, blaKPC, blaNDM, blaOXA, blaCTX-M, vanA, vanB, mcr-1, tetM, ermB, qnrS, aac(6')-Ib-cr
- **42 sequencias FASTA** no NCBI
- **Guias CRISPR** em TSV (colunas: rank, gene, spacer_seq, pam_seq, etc.)
- **Conservation column**: `match` (nao `match_type`)
- **CARD database**: variantes clinicas e antibiotic resistance ontology

---

## Resultados por Experimento

| Task | Descricao | Status | Output Principal |
|------|-----------|--------|------------------|
| A | Guide redesign | PASS | Guias CRISPR otimizados (thermodynamics + PAM-proximal) |
| B | Embedding comparison | PASS | Heatmap Evo2 distances |
| C | Dashboard interativo | PASS | HTML com conservation charts |
| D | Cross-reactivity | PASS | Off-target spacer analysis |
| E | Publication stats | PASS | LaTeX tables + confidence intervals |
| F | Codon usage (RSCU) | PASS | `rscu_analysis.csv` (72KB) |
| G | Hamming distances | PASS | `hamming_distance_matrix.csv` (13KB) |
| H | Coverage heatmap | PASS | HTML gene-family vs match-type |
| I | GC/Tm analysis | PASS | Propriedades biofisicas dos guias |
| J | CARD enrichment | PASS | `card_enrichment.html` (111KB) |
| K | Comprehensive report | PASS | Dashboard multi-secao |
| L | Multiplex optimization | PASS | `multiplex_panel.csv` |
| M | Network adjacency | PASS | `network_adjacency.csv` |
| N | Diversity stats | PASS | `diversity_stats.csv` |
| O | Summary dashboard | PASS | KPI cards + mini-charts |

---

## Principais Outputs

### Na DGX (`~/jv-teste/harvard_bacend/results/`)

| Arquivo | Tamanho | Descricao |
|---------|---------|-----------|
| `card_enrichment.html` | 111KB | Dashboard de variantes CARD enriquecidas |
| `rscu_analysis.csv` | 72KB | Relative Synonymous Codon Usage por familia |
| `hamming_distance_matrix.csv` | 13KB | Distancias entre variantes |
| `multiplex_panel.csv` | — | Painel CRISPR otimizado para multiplex |
| `network_adjacency.csv` | — | Adjacencia gene-gene via similaridade |
| `diversity_stats.csv` | — | Shannon entropy, nucleotide diversity |
| `hamming_distance_matrix.html` | — | Visualizacao interativa |
| `coverage_heatmap.html` | — | Cobertura dos guias por familia |
| `summary_dashboard.html` | — | Overview executivo |

### Pipeline Base (ja no Mac em `/Users/jv/Documents/GitHub/bac_end_harvard/`)

| Arquivo | Descricao |
|---------|-----------|
| `reports/card_integration_report.txt` | Relatorio integracao CARD |
| `reports/card_new_variants.csv` | Novas variantes detectadas |
| `reports/conservation_analysis.tsv` | Analise de conservacao dos guias |
| `reports/covariance_probes/` | 18 features biofisicas |
| `reports/evo2_scoring/` | Scores Evo2 funcionais |
| `reports/oligo_order.tsv` | Ordem de oligos para sintese |
| `reports/panel_report.txt` | Relatorio do painel multiplex |
| `reports/specificity_report.tsv` | Analise de especificidade |
| `reports/tracking_status.csv` | Status de tracking dos genes |

---

## Fixes Aplicados Durante Autoresearch

Durante a execucao inicial, algumas tasks falharam. Os principais problemas e correcoes:

### Task B & K — "No module named 'pandas'"
- **Causa**: Scripts rodaram fora do venv
- **Fix**: `cd /home/oftalmousp/jv-teste/harvard_bacend && source venv/bin/activate`
- Task B ainda tinha "file not found" por caminhos relativos → fix com cd para diretorio correto
- Task K tinha erro de seaborn plot mas exit 0

### Colunas inconsistentes
- `match` foi confundida com `match_type` — ajustado no auto-patcher

### Bug legacy em `card_integration.py` (Mac)
- CSV com chaves inconsistentes entre linhas
- **Fix aplicado no codigo local**: usar `DictWriter` com fieldnames explicitos

---

## Como Baixar os Arquivos da DGX

Quando SSH via cloudflared estiver OK:

```bash
# Copiar todos os outputs do autoresearch
scp -P 2222 oftalmousp@dgx.retina.ia.br:~/jv-teste/harvard_bacend/results/*.csv \
    /Users/jv/Documents/GitHub/bac_end_harvard/results_dgx/

scp -P 2222 oftalmousp@dgx.retina.ia.br:~/jv-teste/harvard_bacend/results/*.html \
    /Users/jv/Documents/GitHub/bac_end_harvard/results_dgx/
```

---

## Proximos Passos

- [ ] Copiar todos os outputs HTML/CSV da DGX quando tunnel estiver OK
- [ ] Incluir `multiplex_panel.csv` no paper/submissao
- [ ] Validar `rscu_analysis` com referencias de literatura
- [ ] Rodar 50 experimentos novos sugeridos (ver `Novos_Experimentos_Gemma4_DGX.md`)

---

## Links

- GitHub local: `/Users/jv/Documents/GitHub/bac_end_harvard/`
- DGX: `~/jv-teste/harvard_bacend/` (acesso via `dgx.retina.ia.br:2222`)
- [[Resumo_Experimentos_DGX_Abril2026]] — Resumo dos 3 experimentos DGX
- [[BlinkTracking_Status]] — Projeto PhD principal

---

*Gerado em 16 de abril de 2026*
*Responsavel: Joao Victor Pacheco Dias*
