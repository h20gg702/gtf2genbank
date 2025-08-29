# gtf2genbank

このリポジトリには、GTF アノテーションとリファレンスゲノム FASTA を **SnapGene Viewer (無料版)** で開ける GenBank (`.gbk`) ファイルに変換するための、個人的なユーティリティを収録しています。  
主に自分の記録・可視化のために作成しましたが、同じようなニーズを持つ方の参考になるかもしれないので公開します。

これらのスクリプトは **本番環境向けに最適化されていません**。  
機能は限定的であり、**利用は自己責任でお願いします！ :)**

## 利用目的

- 複数の転写産物（アイソフォーム）の構造（exon、UTR など）を **素早く可視化**  
- 遺伝子ごとの **転写産物の違いを目で比較**  
- **プライマー設計の補助**（exon 境界やUTRの位置を明確に確認）  

## 使い方

```bash
# 遺伝子シンボルを指定
./gtf2genbank.sh -g TP53 -a gencode.v44.annotation.gtf -f GRCh38.fa

# ENSG ID を指定（バージョン番号は無視）
./gtf2genbank.sh --ensg ENSG00000141510 -a gencode.v44.annotation.gtf -f GRCh38.fa

# ENST ID を指定（バージョン番号は無視）
./gtf2genbank.sh --enst ENST00000413465 -a gencode.v44.annotation.gtf -f GRCh38.fa

# UTR 情報も含める
./gtf2genbank.sh -g TP53 -a gencode.v44.annotation.gtf -f GRCh38.fa -u
```

出力ファイル（`./snapgene_<LABEL>/` に保存）:
- `<LABEL>.gbk` （SnapGene Viewer で開く）
- `<LABEL>.locus.fa`
- `<LABEL>.snapgene.gff3`

## ライセンス

MIT License。本ツールは無保証で提供されます。自己責任でご利用ください。
