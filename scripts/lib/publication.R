`%||%` <- function(a,b) if (is.null(a) || !length(a)) b else a
need <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(absent)) stop('Missing packages: ',paste(absent,collapse=', '),
    '. Run install_figure_dependencies.R, then render again.',call.=FALSE)
}
clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub('_assembly_contigs\\.fasta$', '', x)
  x <- sub('_S[0-9]+_L[0-9]+$', '', x)
  gsub('#','_',x,fixed=TRUE)
}
blank <- function(x) is.na(x) | trimws(as.character(x)) == ''
as_other <- function(x) {
  z<-trimws(as.character(x));z[blank(z)|tolower(z)%in%c('unknown','not known','unavailable')]<-'Other';z
}
stopif <- function(test, text) { if (isTRUE(test)) stop(text,call.=FALSE) }
require_columns <- function(x,cols,what) {
  absent <- setdiff(cols,names(x))
  if(length(absent)) stop(what,': missing columns ',paste(absent,collapse=', '),call.=FALSE)
}
unique_ids <- function(x,what) {
  stopif(any(blank(x)),paste(what,'has empty identifiers'))
  stopif(anyDuplicated(x)>0,paste(what,'has duplicate/normalisation-colliding IDs:',
    paste(head(unique(x[duplicated(x)]),8),collapse=', ')))
  invisible(x)
}
column <- function(x,choices,required=FALSE,default=NA_character_) {
  hit <- choices[choices %in% names(x)]
  if(length(hit)) return(as.character(x[[hit[1]]]))
  if(required) stop('Missing column: expected ',paste(choices,collapse=' / '),call.=FALSE)
  rep(default,nrow(x))
}
numeric_checked <- function(x,what) {
  z <- suppressWarnings(as.numeric(x))
  stopif(any(!blank(x)&is.na(z)),paste('Non-numeric values in',what))
  z
}
read_table <- function(file,header=TRUE) {
  stopif(!file.exists(file),paste('Missing input:',file))
  sep <- if(grepl('\\.tsv$',file,ignore.case=TRUE)) '\t' else ','
  x <- utils::read.table(file,header=header,sep=sep,quote='"',comment.char='',
    stringsAsFactors=FALSE,colClasses='character',check.names=FALSE,
    na.strings=c('','NA'),fill=FALSE)
  stopif(anyDuplicated(names(x))>0,paste('Duplicate column names in',file))
  x
}
identify <- function(x,choices=c('isolate_id','Isolate_ID','Accession','Isolate','Isolates')) {
  x$original_isolate_id <- column(x,choices,TRUE)
  x$isolate_id <- clean_id(x$original_isolate_id)
  unique_ids(x$isolate_id,'Metadata'); x
}
coalesce_checked <- function(a,b,what) {
  known <- !blank(a)&!blank(b)
  stopif(any(known&as.character(a)!=as.character(b)),paste('Conflicting',what,
    '- reconcile the source tables instead of silently selecting one.'))
  a[blank(a)] <- b[blank(a)]; a
}
merge_checked <- function(a,b,columns,what) {
  unique_ids(b$isolate_id,what); ii<-match(a$isolate_id,b$isolate_id)
  for(n in intersect(columns,names(b))) {
    v<-b[[n]][ii]
    if(n %in% names(a)) a[[n]]<-coalesce_checked(a[[n]],v,paste(what,n)) else a[[n]]<-v
  };a
}
validate_distance <- function(D) {
  stopif(!is.matrix(D)||!is.numeric(D)||nrow(D)!=ncol(D),'SNP matrix must be numeric and square')
  stopif(is.null(rownames(D))||is.null(colnames(D)),'SNP matrix must name both axes')
  rownames(D)<-clean_id(rownames(D));colnames(D)<-clean_id(colnames(D))
  unique_ids(rownames(D),'SNP rows');unique_ids(colnames(D),'SNP columns')
  stopif(!setequal(rownames(D),colnames(D)),'SNP row/column sets differ')
  D<-D[,match(rownames(D),colnames(D)),drop=FALSE]
  stopif(any(!is.finite(D)),'Missing/non-numeric SNP values are not zero')
  stopif(any(D<0)||any(abs(D-round(D))>1e-8),'SNP counts must be non-negative integers')
  stopif(any(diag(D)!=0)||max(abs(D-t(D)))>1e-8,'SNP matrix is not symmetric with a zero diagonal')
  D
}
parse_date_checked <- function(x) {
  x<-trimws(as.character(x));ans<-as.Date(rep(NA_character_,length(x)))
  for(fmt in c('%Y-%m-%d','%d/%m/%Y','%d-%m-%Y','%d %b %Y','%d %B %Y')) {
    ii<-which(is.na(ans)&!blank(x));if(!length(ii))break

    candidate<-suppressWarnings(as.Date(x[ii],format=fmt))
    good<-!is.na(candidate)&format(candidate,fmt)==x[ii]
    ans[ii[good]]<-candidate[good]
  }; ans
}

INK <- '#21313F'; MUTED <- '#647784'; GRID <- '#DCE4E8'; PAPER <- '#FFFFFF'
PAL <- c('Malawi'='#AF493C','South Africa'='#2C6D91','Mozambique'='#B58A32',
  'Nigeria'='#5B8354','United Kingdom'='#8064A2','Brazil'='#BA7752',
  'United States'='#6989A6','Taiwan'='#A66888','Mexico'='#548F83',
  'Blood culture'='#AF493C','CSF'='#8064A2','River water'='#3186A0',
  'Neonatal-unit sampling'='#56856B','Other'='#C2CAD0')
colours <- function(values) {
  values<-sort(unique(as_other(values)))
  extra<-setdiff(values,names(PAL))
  more<-setNames(vapply(extra,function(z){
    h<-sum(as.integer(charToRaw(enc2utf8(z)))*seq_along(charToRaw(enc2utf8(z))))%%360
    grDevices::hcl(h,50,52)
  },character(1)),extra)
  c(PAL,more)[values]
}
pub_theme <- function(base=9) ggplot2::theme_minimal(base_size=base,base_family='sans')+
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
    panel.grid.major.x=ggplot2::element_blank(),axis.text=ggplot2::element_text(colour=INK),
    axis.title=ggplot2::element_text(colour=INK),
    plot.caption=ggplot2::element_text(size=7,colour=MUTED,hjust=0),
    legend.position='bottom',legend.title=ggplot2::element_text(face='bold',size=8),
    legend.text=ggplot2::element_text(size=7),
    legend.key.height=grid::unit(3,'mm'),
    plot.margin=ggplot2::margin(6,7,6,7),
    strip.text=ggplot2::element_text(face='bold',colour=INK))

native <- function(w,h) grid::pushViewport(grid::viewport(xscale=c(0,w),yscale=c(0,h)))
gtxt <- function(label,x,y,size=8,col=INK,just='left',rot=0,bold=FALSE) {
  grid::grid.text(label,x=x,y=y,default.units='native',just=just,rot=rot,
    gp=grid::gpar(fontfamily='sans',fontsize=size,col=col,fontface=if(bold)'bold' else 'plain'))
}
gline <- function(x0,y0,x1,y1,col=GRID,lwd=.5) grid::grid.segments(x0,y0,x1,y1,
  default.units='native',gp=grid::gpar(col=col,lwd=lwd,lineend='round'))
grect <- function(x,y,w,h,fill=NA,col=NA,lwd=.3) grid::grid.rect(x,y,w,h,
  default.units='native',gp=grid::gpar(fill=fill,col=col,lwd=lwd))
page_locator <- function(w,h,page=NULL) {

  if(!is.null(page))gtxt(page,w-9,h-6,size=7,just='right',col=MUTED)
  invisible(NULL)
}
wrapped <- function(x,width=130) paste(strwrap(x,width=width),collapse='\n')
wrap_label <- function(x,width_mm,size=7) {

  vapply(x,function(s){
    if(grid::convertWidth(grid::grobWidth(grid::textGrob(s,gp=grid::gpar(fontsize=size))),
       'mm',valueOnly=TRUE)<=width_mm)return(s)
    paste(strwrap(s,width=max(12,floor(width_mm/(size*.18)))),collapse='\n')
  },character(1))
}

default_paths <- function() c(
  master='data/data_clean/all_good_isangi_genomes_metadata_amr_ebg.csv',
  genomes='data/data_raw/good_genomes.csv',metadata='data/data_raw/metadata.csv',
  tree='data/data_raw/RAxML_best_trees/RAxML_bestTree.ebg25_quality_assured_tree',
  snp_matrix='data/data_clean/ebg_25_snpmatrix.tsv',
  cc25_meta='data/data_clean/EBG_25_metadata.csv',
  mlst='data/data_clean/all_isolates_mlst_24052024.csv',
  lineages='data/data_clean/ebg25.snp-sites.fastbaps_lineages.csv',
  classes='data/data_clean/all_isangi_assemblies_amrfinder2_20062024.number_of_classes.csv',
  malawi='data/data_clean/malawi_metadata_for_plotting.csv',
  amr='publication_inputs/amr_determinants.csv',
  replicons='publication_inputs/plasmid_replicons.csv',
  aliases='publication_inputs/accession_aliases.csv',
  annotations='publication_inputs/isolate_annotations.csv',
  map='publication_inputs/map_layers.gpkg',
  site_members='publication_inputs/site_isolate_membership.csv',
  burden='publication_inputs/ncbi_amr_counts.csv',
  plasmids='publication_inputs/plasmids.csv',
  features='publication_inputs/plasmid_features.csv',
  regions='publication_inputs/plasmid_regions.csv',
  links='publication_inputs/plasmid_alignments.csv')
init <- function(params,bundle_dir=getwd()) {
  s<-new.env(parent=emptyenv());s$cfg<-params
  s$root<-normalizePath(params$project_root,mustWork=TRUE)
  s$bundle<-normalizePath(bundle_dir,mustWork=TRUE);s$paths<-default_paths()
  override<-params$paths_file %||% ''
  if(nzchar(override)) {
    p<-if(grepl('^(/|[A-Za-z]:)',override))override else file.path(s$root,override)
    o<-read_table(p);require_columns(o,c('key','path'),'Path overrides');unique_ids(o$key,'Path keys')
    stopif(any(!o$key%in%names(s$paths)),'Unknown path override key')
    s$paths[o$key]<-o$path
  }
  stamp<-paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_',Sys.getpid())
  s$out<-file.path(s$root,'outputs','publication_figures',stamp)
  folders<-c('main','supplementary','source_data','audit','captions','_failed')
  if(isTRUE(params$vector_archive))folders<-c(folders,'vector_archive/main','vector_archive/supplementary')
  for(d in folders)
    dir.create(file.path(s$out,d),recursive=TRUE,showWarnings=FALSE)
  s$status<-data.frame(panel=character(),state=character(),detail=character())
  s$warnings<-character();s$registry<-NULL;s$cc<-NULL;s$master<-NULL;s$snp<-NULL
  s$used<-character()
  manifest<-data.frame(key=names(s$paths),path=unname(s$paths),
    present=vapply(names(s$paths),function(k)file.exists(input_path(s,k)),logical(1)))
  utils::write.csv(manifest,file.path(s$out,'audit','input_manifest.csv'),row.names=FALSE)
  writeLines(capture.output(dput(params)),file.path(s$out,'audit','parameters.R'))
  writeLines(c('Revision: 2026-09-19-final-paper-ready',
    'Main 600-dpi TIFFs: main/',
    'Supplementary 600-dpi TIFFs: supplementary/',
    'Numeric matrices, isolate metadata, memberships and insertion guides: source_data/',
    'The terminal prints this timestamped output directory when rendering finishes.'),
    file.path(s$out,'WHERE_ARE_MY_FIGURES.txt'))
  writeLines('2026-09-19-final-paper-ready',file.path(s$out,'audit','FIGURE_BUILD_REVISION.txt'))
  s
}
input_path <- function(s,key) {
  p<-s$paths[[key]]
  if(grepl('^(/|[A-Za-z]:)',p))p else file.path(s$root,p)
}
has <- function(s,key) file.exists(input_path(s,key))
read_input <- function(s,key,header=TRUE) {
  p<-input_path(s,key);x<-read_table(p,header);s$used<-unique(c(s$used,p));x
}
read_bundled <- function(s,name) {
  p<-file.path(s$bundle,'figure_inputs',name)
  x<-read_table(p);s$used<-unique(c(s$used,p));x
}
export <- function(s,x,name) {
  utils::write.csv(x,file.path(s$out,'source_data',paste0(name,'.csv')),row.names=FALSE,na='')
}
caption <- function(s,name,text) writeLines(text,file.path(s$out,'captions',paste0(name,'.txt')))
note <- function(s,text) {s$warnings<-unique(c(s$warnings,text));message(text)}
save_pdf_drawing <- function(file,draw,w,h) {
  grDevices::pdf(file,width=w/25.4,height=h/25.4,useDingbats=FALSE)
  tryCatch({grid::grid.newpage();draw(w,h)},finally=grDevices::dev.off())
}
save_tiff_drawing <- function(s,file,draw,w,h) {
  need('ragg')
  dpi<-as.integer(s$cfg$tiff_dpi)
  stopif(!is.finite(dpi)||dpi<300L,'tiff_dpi must be at least 300')
  ragg::agg_tiff(file,width=w,height=h,units='mm',res=dpi,compression='lzw',
    background='white',bitsize=8)
  tryCatch({grid::grid.newpage();draw(w,h)},finally=grDevices::dev.off())
}
save_drawing <- function(s,name,draw,w,h,supp=FALSE,pdf=TRUE) {
  folder<-if(supp)'supplementary' else 'main';base<-file.path(s$out,folder,name)
  if(isTRUE(s$cfg$vector_archive)&&isTRUE(pdf))
    save_pdf_drawing(file.path(s$out,'vector_archive',folder,paste0(name,'.pdf')),draw,w,h)
  save_tiff_drawing(s,paste0(base,'.tiff'),draw,w,h)
  invisible(base)
}
save_gg <- function(s,p,name,w=180,h=140,supp=FALSE) {
  save_drawing(s,name,function(w,h)print(p,newpage=FALSE),w,h,supp)
}
run_panel <- function(s,name,fun) {
  message('\n',name)
  s$warnings<-character();before<-list.files(s$out,recursive=TRUE,full.names=TRUE)
  state<-'BUILT';detail<-''
  tryCatch(withCallingHandlers(fun(),warning=function(w){
    s$warnings<-unique(c(s$warnings,conditionMessage(w)));invokeRestart('muffleWarning')
  }),error=function(e){state<<-'BLOCKED';detail<<-conditionMessage(e)})
  if(state=='BUILT'&&length(s$warnings)) {state<-'BUILT - REVIEW';detail<-paste(s$warnings,collapse=' | ')}
  if(state=='BLOCKED') {

    after<-setdiff(list.files(s$out,recursive=TRUE,full.names=TRUE),before)
    graphics<-after[grepl('\\.(pdf|tiff)$',after,ignore.case=TRUE)]
    for(f in graphics)file.rename(f,file.path(s$out,'_failed',paste0(gsub('[^A-Za-z0-9]','_',name),'_',basename(f))))
  }
  s$status<-rbind(s$status,data.frame(panel=name,state=state,detail=detail))
  utils::write.csv(s$status,file.path(s$out,'audit','build_status.csv'),row.names=FALSE)
  message(state,if(nzchar(detail))paste0(': ',detail) else '')
}

valid_run <- function(x) !blank(x)&grepl('^[ESD]RR[0-9]+(;[ESD]RR[0-9]+)*$',x)
canonical_runs <- function(x) vapply(x,function(z){
  if(blank(z))return(NA_character_)
  z<-gsub('[[:space:]]','',z);paste(sort(unique(strsplit(z,';',fixed=TRUE)[[1]])),collapse=';')
},character(1),USE.NAMES=FALSE)
load_registry <- function(s) {
  r<-read_bundled(s,'ena_accession_crosswalk.csv');require_columns(r,c('isolate_id','run_accession'),'ENA register')
  r$isolate_id<-clean_id(r$isolate_id);r$run_accession<-canonical_runs(r$run_accession)
  unique_ids(r$isolate_id,'ENA register');stopif(any(!valid_run(r$run_accession)),'Invalid run accession in register')
  if(has(s,'aliases')) {
    a<-identify(read_input(s,'aliases'));require_columns(a,'run_accession','Additional accessions')
    a$run_accession<-canonical_runs(a$run_accession)
    stopif(any(!valid_run(a$run_accession)),'Aliases require real run accessions')
    for(i in seq_len(nrow(a))) {
      k<-match(a$isolate_id[i],r$isolate_id)
      if(!is.na(k))r$run_accession[k]<-coalesce_checked(r$run_accession[k],a$run_accession[i],'accession mapping')
      else {
        new<-as.data.frame(setNames(rep(list(NA_character_),ncol(r)),names(r)),stringsAsFactors=FALSE)
        for(j in intersect(names(r),names(a)))new[[j]]<-a[[j]][i]
        r<-rbind(r,new)
      }
    }
  }

  rr<-unlist(strsplit(r$run_accession,';',fixed=TRUE),use.names=FALSE)
  unique_ids(rr,'Run-to-isolate map')
  s$registry<-r;export(s,r,'accession_register_80_plus_curated_aliases');r
}
attach_accessions <- function(s,m) {
  if(is.null(s$registry))load_registry(s)
  unique_ids(m$isolate_id,'Plot metadata')
  ii<-match(m$isolate_id,s$registry$isolate_id)
  given<-canonical_runs(column(m,c('run_accession','ena_run','Run','SRA_Accession','SRA_accession','SRA')))
  bad<-!blank(given)&!valid_run(given)
  if(any(bad)) {note(s,'Some metadata accession fields are not read runs; they were retained in source data, not used as run labels.');given[bad]<-NA_character_}
  inherent<-ifelse(valid_run(m$isolate_id),m$isolate_id,NA_character_)
  m$run_accession<-coalesce_checked(coalesce_checked(given,s$registry$run_accession[ii],'ENA/run mapping'),inherent,'intrinsic run IDs')

  same_as_run<-!blank(m$run_accession)&m$isolate_id==m$run_accession
  m$plot_label<-ifelse(same_as_run,m$isolate_id,
    paste(m$isolate_id,ifelse(blank(m$run_accession),'RUN UNRESOLVED',m$run_accession),sep=' | '))
  missing<-blank(m$run_accession)
  if(any(missing)) {
    export(s,m[missing,c('isolate_id','plot_label'),drop=FALSE],paste0('unresolved_runs_',nrow(m)))
    text<-paste(sum(missing),'isolates have unresolved read-run accessions.')
    if(isTRUE(s$cfg$strict))stop(text,call.=FALSE) else note(s,text)
  };m
}
canonical_meta <- function(x) {
  z<-identify(x)
  z$sequence_type<-column(z,c('sequence_type','ST','V3','MLST'))
  z$country<-column(z,c('country','Country'))
  z$source<-column(z,c('source','source_class','Comment','Comment.x','Class'))
  z$year<-column(z,c('year','Collection Year','Collection.Year','Collection.Year.x'))
  z$lineage<-column(z,c('lineage','Level 1','Level.1'))
  z$amr_count<-column(z,c('amr_count','Unique_Class_Count'))
  z
}
load_master <- function(s) {
  if(!is.null(s$master))return(s$master)
  if(has(s,'master')) m<-canonical_meta(read_input(s,'master')) else {
    stopif(!has(s,'genomes')||!has(s,'metadata'),'No master metadata: provide cleaned master or good_genomes.csv + metadata.csv')
    g<-canonical_meta(read_input(s,'genomes'));d<-canonical_meta(read_input(s,'metadata'))
    m<-merge_checked(g,d,setdiff(names(d),c('isolate_id','original_isolate_id')),'raw metadata')
  }
  m<-attach_accessions(s,m);s$master<-m;export(s,m,'all_genomes_metadata');m
}
load_tree_set <- function(s) {
  if(!is.null(s$cc))return(s$cc)
  need(c('ape','phangorn'))
  stopif(!has(s,'tree'),paste('Missing CC25 Newick:',s$paths[['tree']]))
  tr<-ape::read.tree(input_path(s,'tree'));s$used<-unique(c(s$used,input_path(s,'tree')))
  stopif(!inherits(tr,'phylo'),'Expected exactly one phylogenetic tree')
  original_tip<-tr$tip.label;tr$tip.label<-clean_id(tr$tip.label);unique_ids(tr$tip.label,'Tree tips')
  stopif(is.null(tr$edge.length)||any(!is.finite(tr$edge.length))||any(tr$edge.length<0),
    'Tree must have finite non-negative branch lengths')
  before<-ape::cophenetic.phylo(tr)
  if(isTRUE(s$cfg$midpoint_root))tr<-phangorn::midpoint(tr)
  tr<-ape::ladderize(tr,right=TRUE)
  after<-ape::cophenetic.phylo(tr)[rownames(before),colnames(before)]
  stopif(max(abs(before-after))>1e-8*max(1,max(before)),'Rerooting changed pairwise patristic distances')
  m<-canonical_meta(read_input(s,'cc25_meta'))
  stopif(!all(tr$tip.label%in%m$isolate_id),'Some tree tips have no EBG_25 metadata row')
  extra<-setdiff(m$isolate_id,tr$tip.label)
  export(s,data.frame(isolate_id=extra),'CC25_metadata_not_in_tree')
  m<-m[match(tr$tip.label,m$isolate_id),,drop=FALSE]
  if(has(s,'mlst')) {
    raw<-read_input(s,'mlst',FALSE);stopif(ncol(raw)<3,'MLST file requires at least three columns')
    ml<-data.frame(isolate_id=clean_id(raw[[1]]),sequence_type=raw[[3]])
    m<-merge_checked(m,ml,'sequence_type','MLST')
  }
  if(has(s,'lineages')) {
    raw<-identify(read_input(s,'lineages'));lc<-data.frame(isolate_id=raw$isolate_id,lineage=column(raw,c('Level 1','Level.1','lineage'),TRUE))
    m<-merge_checked(m,lc,'lineage','hierBAPS')
  }
  if(has(s,'classes')) {
    raw<-identify(read_input(s,'classes'));ac<-data.frame(isolate_id=raw$isolate_id,amr_count=column(raw,c('Unique_Class_Count','amr_count'),TRUE))
    m<-merge_checked(m,ac,'amr_count','AMR counts')
  }
  if(has(s,'annotations')) {
    a<-canonical_meta(read_input(s,'annotations'))
    m<-merge_checked(m,a,c('sequence_type','country','source','year','lineage','amr_count'),'curated annotations')
  }
  m$amr_count<-numeric_checked(m$amr_count,'AMR subclass counts')
  stopif(any(!is.na(m$amr_count)&(m$amr_count<0|m$amr_count!=round(m$amr_count))),'AMR subclass counts must be non-negative integers')
  m<-attach_accessions(s,m)
  m$original_tree_tip<-original_tip[match(m$isolate_id,clean_id(original_tip))]
  s$cc<-list(tree=tr,meta=m)
  export(s,m,'CC25_exact_tip_metadata')
  if(length(tr$tip.label)!=327L)note(s,paste('CC25 has',length(tr$tip.label),'tips; the manuscript reports 327. No tips were forced or discarded.'))
  if(any(blank(m$sequence_type)))note(s,'Some sequence types are unavailable; they remain on CC25 under Other and cannot be assigned to the ST335 subset.')
  s$cc
}

load_snp_matrix <- function(s) {
  if(!is.null(s$snp))return(s$snp)
  stopif(!has(s,'snp_matrix'),paste('Missing pairwise SNP matrix:',s$paths[['snp_matrix']]))
  p<-input_path(s,'snp_matrix')
  raw<-utils::read.delim(p,header=TRUE,check.names=FALSE,row.names=NULL,quote='',comment.char='',
    stringsAsFactors=FALSE,colClasses='character');s$used<-unique(c(s$used,p))
  stopif(ncol(raw)!=nrow(raw)+1L,
    'SNP TSV must contain one isolate-ID field followed by one named numeric column per isolate.')
  ids<-raw[[1L]];values<-raw[-1L]
  ids<-clean_id(ids);unique_ids(ids,'SNP matrix rows')
  names(values)<-clean_id(names(values));unique_ids(names(values),'SNP matrix columns')
  D<-as.matrix(data.frame(lapply(values,numeric_checked,what='pairwise SNP counts'),check.names=FALSE))
  rownames(D)<-ids;colnames(D)<-names(values);D<-validate_distance(D)
  s$snp<-D;D
}
export_snp_matrix <- function(s,D,prefix) {
  wide<-data.frame(isolate_id=rownames(D),as.data.frame(D,check.names=FALSE),
    check.names=FALSE,stringsAsFactors=FALSE)
  export(s,wide,paste0(prefix,'_pairwise_SNP_matrix_tree_order'))
  utils::write.table(wide,file.path(s$out,'source_data',paste0(prefix,'_pairwise_SNP_matrix_tree_order.tsv')),
    sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE,na='')
  pair<-which(upper.tri(D),arr.ind=TRUE)
  long<-data.frame(isolate_a=rownames(D)[pair[,1]],isolate_b=colnames(D)[pair[,2]],
    snp_distance=as.integer(D[pair]),stringsAsFactors=FALSE)
  export(s,long,paste0(prefix,'_pairwise_SNP_long'))
  nonzero<-long$snp_distance[long$snp_distance>0]
  summary<-data.frame(matrix=prefix,n_isolates=nrow(D),n_pairs=nrow(long),
    minimum=if(length(long$snp_distance))min(long$snp_distance) else NA_real_,
    minimum_nonzero=if(length(nonzero))min(nonzero) else NA_real_,
    median=if(length(long$snp_distance))stats::median(long$snp_distance) else NA_real_,
    maximum=if(length(long$snp_distance))max(long$snp_distance) else NA_real_)
  export(s,summary,paste0(prefix,'_pairwise_SNP_summary'))
  invisible(D)
}
safe_slug <- function(x) {
  z<-gsub('[^A-Za-z0-9]+','_',trimws(as.character(x)));z<-gsub('^_+|_+$','',z)
  ifelse(nzchar(z),z,'group')
}
matrix_track_bundle <- function(meta,scope,status=NULL) {
  unique_ids(meta$isolate_id,paste(scope,'SNP metadata IDs'))
  tk<-track_data(meta)
  keep<-if(scope=='ST335')c('Lineage','Country','Source','Year') else
    c('ST','Lineage','Country','Source','Year','AMR subclass count')
  keep<-intersect(keep,names(tk$data));tk$data<-tk$data[keep];tk$cols<-tk$cols[keep]
  export_data<-tk$data
  if(!is.null(status)) {
    components<-attr(status,'component_status')
    unique_ids(status$isolate_id,paste(scope,'SNP status IDs'))
    status<-status[match(meta$isolate_id,status$isolate_id),,drop=FALSE]
    stopif(anyNA(status$isolate_id),paste(scope,'SNP metadata lacks status rows'))
    shown<-intersect(c('XDR status','CARB-R / MAC-R'),names(status))
    for(k in shown) {
      tk$data[[k]]<-as_other(status[[k]])
      tk$cols[[k]]<-if(k=='XDR status')
        c(`non-XDR`='#D7D9DB',XDR='#F05A61',Other='#FFFFFF') else
        c(`non-CARB-R / MAC-R`='#D7D9DB',`CARB-R / MAC-R`='#F05A61',Other='#FFFFFF')
      tk$cols[[k]]<-tk$cols[[k]][names(tk$cols[[k]])%in%unique(tk$data[[k]])]
      export_data[[k]]<-tk$data[[k]]
    }
    if(!is.null(components)) {
      components<-components[match(meta$isolate_id,components$isolate_id),,drop=FALSE]
      stopif(anyNA(components$isolate_id),paste(scope,'SNP metadata lacks component-status rows'))
      for(k in setdiff(names(components),'isolate_id'))export_data[[k]]<-as_other(components[[k]])
    }
  }
  list(data=tk$data,cols=tk$cols,export_data=export_data)
}
snp_order_metadata <- function(meta,tk) {
  data.frame(matrix_order=seq_len(nrow(meta)),isolate_id=meta$isolate_id,
    run_accession=if('run_accession'%in%names(meta))meta$run_accession else NA_character_,
    plot_label=if('plot_label'%in%names(meta))meta$plot_label else meta$isolate_id,
    recorded_sequence_type=as_other(meta$sequence_type),
    recorded_lineage=as_other(meta$lineage),recorded_country=as_other(meta$country),
    recorded_source=as_other(meta$source),recorded_year=as_other(meta$year),
    recorded_amr_subclass_count=as_other(meta$amr_count),
    tk$export_data,check.names=FALSE,stringsAsFactors=FALSE)
}
snp_heatmap_draw <- function(D,tk,label_ids=FALSE,common_max=NULL) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    stopif(!identical(rownames(D),colnames(D)),'SNP heatmap axes are not in identical order')
    n<-nrow(D);ntracks<-ncol(tk$data);trackw<-2.15;legendw<-76
    label_left<-if(label_ids)45 else 9;bottom<-if(label_ids)46 else 10
    x0<-label_left+ntracks*trackw
    side<-min(w-legendw-x0-7,h-bottom-ntracks*trackw-7)
    stopif(side<70,'SNP heatmap canvas is too small for the metadata and legend panels')
    y0<-bottom;cell<-side/n
    high<-common_max%||%max(D);if(!is.finite(high)||high<0)high<-0
    pal<-grDevices::colorRampPalette(c('#F7FAFC','#D5E7EA','#86B9C2','#347F91','#163F5A'))(256)
    scaled<-if(high==0)matrix(0,nrow(D),ncol(D)) else sqrt(D/high)
    idx<-pmax(1L,pmin(256L,round(scaled*255)+1L))
    fill<-matrix(pal[idx],nrow=nrow(D),ncol=ncol(D))
    grid::grid.raster(fill,x=x0+side/2,y=y0+side/2,width=side,height=side,
      default.units='native',interpolate=FALSE)
    grect(x0+side/2,y0+side/2,side,side,fill=NA,col='#596B76',lwd=.5)
    xc<-x0+(seq_len(n)-.5)*cell;yc<-y0+side-(seq_len(n)-.5)*cell
    short<-c(ST='ST',Lineage='Lineage',Country='Country',Source='Source',Year='Year',
      `AMR subclass count`='AMR count',`XDR status`='XDR',`CARB-R / MAC-R`='CARB+MAC')
    for(j in seq_len(ntracks)) {
      k<-names(tk$data)[j];v<-tk$data[[k]];cols<-tk$cols[[k]];tile<-unname(cols[v])
      grect(xc,y0+side+(j-.5)*trackw,cell*1.02,trackw*.88,fill=tile,col=NA)
      grect(x0-(j-.5)*trackw,yc,trackw*.88,cell*1.02,fill=tile,col=NA)
      lab<-if(k%in%names(short))unname(short[[k]]) else k
      gtxt(lab,x0-1,y0+side+(j-.5)*trackw,size=5.2,just='right',bold=TRUE)
    }
    if(label_ids) {
      stopif(n>45L,'Labelled SNP audit blocks are limited to 45 isolates for legibility')
      fsize<-if(n<=25)5.4 else 4.7
      gtxt(rownames(D),x0-ntracks*trackw-1,yc,size=fsize,just='right')
      gtxt(colnames(D),xc,y0-2,size=fsize,just='right',rot=58)
    }
    legend_x<-x0+side+8
    legend_bottom<-draw_side_legends(tk,legend_x,h-7,w-legend_x-4,
      heading_size=6.1,item_size=5.2,pitch=2.65,gap=2.6)
    stopif(legend_bottom<25,'SNP metadata legend overflow: reduce displayed categories')
    cbx<-legend_x;cby<-11;cbw<-min(50,w-legend_x-7);cbh<-3.2
    for(i in seq_len(128))grect(cbx+(i-.5)*cbw/128,cby,cbw/128+.03,cbh,
      fill=pal[round(1+(i-1)*255/127)])
    tick<-seq(0,1,length.out=5);raw<-round(tick^2*high)
    for(i in seq_along(tick))gtxt(format(raw[i],big.mark=','),cbx+tick[i]*cbw,cby-3.3,
      size=5.2,just=if(i==1)'left' else if(i==length(tick))'right' else 'centre')
    gtxt('Pairwise SNP distance (square-root colour scale)',cbx,cby+5.0,size=5.6)
  }
}
write_snp_outputs <- function(s,ids,prefix,require_exact=FALSE,meta,status=NULL,
  audit_extras=identical(prefix,'CC25')) {
  D<-load_snp_matrix(s)
  missing<-setdiff(ids,rownames(D));extra<-setdiff(rownames(D),ids)
  stopif(length(missing)>0,paste(prefix,'tree contains isolates absent from the SNP matrix:',
    paste(head(missing,8),collapse=', ')))
  if(length(extra)&&isTRUE(audit_extras)) {
    export(s,data.frame(isolate_id=extra,
      exclusion_reason=paste('Present in supplied SNP matrix but outside',prefix,'plotted tree')),
      paste0(prefix,'_SNP_matrix_isolates_outside_tree'))
    note(s,paste(prefix,'SNP matrix contains',length(extra),
      'additional isolates outside the plotted tree; they are listed in source_data and explicitly excluded from this tree-matched matrix.'))
  }
  if(isTRUE(require_exact))stopif(anyDuplicated(ids)>0,'Tree tip identifiers are duplicated after normalisation')
  D<-D[ids,ids,drop=FALSE];D<-validate_distance(D)
  unique_ids(meta$isolate_id,paste(prefix,'SNP metadata IDs'))
  meta<-meta[match(rownames(D),meta$isolate_id),,drop=FALSE]
  stopif(anyNA(meta$isolate_id)||!identical(meta$isolate_id,rownames(D)),
    paste(prefix,'SNP metadata could not be aligned exactly to matrix order'))
  tk<-matrix_track_bundle(meta,prefix,status)
  order_meta<-snp_order_metadata(meta,tk)
  export(s,order_meta,paste0(prefix,'_SNP_isolate_order_metadata'))
  export_snp_matrix(s,D,prefix)
  save_drawing(s,paste0('Supplement_',prefix,'_pairwise_SNP_heatmap'),
    snp_heatmap_draw(D,tk,FALSE,max(D)),297,210,supp=TRUE)
  invisible(list(D=D,meta=meta,tracks=tk))
}
write_snp_panel <- function(s,parent_D,parent_meta,ids,panel,rule,scope,status=NULL,
  common_max=max(parent_D),label_ids=length(ids)<=45L) {
  ids<-rownames(parent_D)[rownames(parent_D)%in%ids]
  if(length(ids)<2L)return(list(manifest=data.frame(panel=panel,rule=rule,n=length(ids),
    labelled=FALSE,state='SKIPPED: fewer than two verified isolates',filename='',stringsAsFactors=FALSE),
    membership=data.frame()))
  D<-validate_distance(parent_D[ids,ids,drop=FALSE])
  meta<-parent_meta[match(ids,parent_meta$isolate_id),,drop=FALSE]
  stopif(anyNA(meta$isolate_id)||!identical(meta$isolate_id,rownames(D)),
    paste(panel,'metadata failed exact matrix-order alignment'))
  tk<-matrix_track_bundle(meta,scope,status)
  slug<-safe_slug(panel);figure<-paste0('Supplement_SNP_',slug)
  export_snp_matrix(s,D,paste0('SNP_',slug))
  order_meta<-snp_order_metadata(meta,tk)
  export(s,order_meta,paste0('SNP_',slug,'_isolate_order_metadata'))
  save_drawing(s,figure,snp_heatmap_draw(D,tk,label_ids,common_max),260,190,supp=TRUE)
  caption(s,figure,paste('Direct subset of the supplied',scope,
    'pairwise SNP matrix in inherited phylogenetic order.',rule,
    if(label_ids)'Both axes print isolate IDs.' else
      'Axis IDs are omitted at this scale; the exact order and isolate metadata are supplied in source_data.',
    'Top and left strips repeat the same isolate metadata. Other denotes unavailable or uncategorised metadata.',
    'The square-root colour scale is labelled in raw SNP counts and uses the maximum of the corresponding full parent matrix.'))
  list(manifest=data.frame(panel=panel,rule=rule,n=length(ids),labelled=label_ids,
    state='BUILT',filename=paste0(figure,'.tiff'),stringsAsFactors=FALSE),
    membership=data.frame(panel=panel,matrix_order=seq_along(ids),isolate_id=ids,
      inclusion_rule=rule,stringsAsFactors=FALSE))
}
write_snp_labelled_blocks <- function(s,D,meta,prefix,scope,status=NULL,block_size=40L) {
  starts<-seq.int(1,nrow(D),by=block_size);out<-vector('list',length(starts))
  for(i in seq_along(starts)) {
    rows<-seq.int(starts[i],min(nrow(D),starts[i]+block_size-1L))
    out[[i]]<-write_snp_panel(s,D,meta,rownames(D)[rows],
      paste0(prefix,'_labelled_block_',sprintf('%02d',i)),
      sprintf('Tree-order audit block %d: matrix positions %d-%d of %d; this is a labelled principal submatrix in inherited order, not an independently clustered subset, and cross-block distances remain in the authoritative full matrix.',
        i,min(rows),max(rows),nrow(D)),scope,status,max(D),TRUE)
  }
  out
}
pairwise_group_comparison <- function(D,groups,analysis) {
  groups<-as_other(groups);names(groups)<-rownames(D)
  keep<-groups!='Other';D<-D[keep,keep,drop=FALSE];groups<-groups[keep]
  if(nrow(D)<2L)return(list(long=data.frame(),summary=data.frame()))
  pair<-which(upper.tri(D),arr.ind=TRUE);ga<-groups[pair[,1]];gb<-groups[pair[,2]]
  swap<-ga>gb;group_a<-ifelse(swap,gb,ga);group_b<-ifelse(swap,ga,gb)
  id_a<-rownames(D)[pair[,1]];id_b<-colnames(D)[pair[,2]]
  isolate_a<-ifelse(swap,id_b,id_a);isolate_b<-ifelse(swap,id_a,id_b)
  long<-data.frame(analysis=analysis,group_a=group_a,group_b=group_b,
    pair_class=ifelse(group_a==group_b,paste('Within',group_a),
      paste('Between',group_a,'and',group_b)),
    isolate_a=isolate_a,isolate_b=isolate_b,
    snp_distance=as.integer(D[pair]),stringsAsFactors=FALSE)
  split_rows<-split(seq_len(nrow(long)),paste(long$group_a,long$group_b,sep='\t'))
  summary<-do.call(rbind,lapply(split_rows,function(ii){
    z<-long$snp_distance[ii];a<-long$group_a[ii[1]];b<-long$group_b[ii[1]]
    q<-stats::quantile(z,c(.25,.5,.75),names=FALSE,type=7)
    data.frame(analysis=analysis,group_a=a,group_b=b,pair_class=long$pair_class[ii[1]],
      n_isolates_a=sum(groups==a),n_isolates_b=sum(groups==b),n_unique_pairs=length(z),
      minimum=min(z),q1=q[1],median=q[2],q3=q[3],maximum=max(z),
      n_le_0=sum(z<=0),n_le_2=sum(z<=2),n_le_5=sum(z<=5),n_le_10=sum(z<=10),
      stringsAsFactors=FALSE)
  }))
  list(long=long,summary=summary)
}

tree_coordinates <- function(tr) {
  nt<-length(tr$tip.label); nn<-nt+tr$Nnode
  root<-setdiff(tr$edge[,1],tr$edge[,2]);stopif(length(root)!=1,'Tree root is not unique')
  child<-split(tr$edge[,2],tr$edge[,1]);order<-integer();ys<-rep(NA_real_,nn)
  visit<-function(v) {
    if(v<=nt) {order<<-c(order,v);return(invisible(NULL))}
    for(k in child[[as.character(v)]])visit(k)
  };visit(root)
  ys[order]<-seq_along(order)
  yy<-function(v){if(is.na(ys[v]))ys[v]<<-mean(vapply(child[[as.character(v)]],yy,numeric(1)));ys[v]}
  yy(root);xs<-ape::node.depth.edgelength(tr)
  horizontal<-data.frame(x0=xs[tr$edge[,1]],x1=xs[tr$edge[,2]],r0=ys[tr$edge[,2]],r1=ys[tr$edge[,2]])
  vertical<-do.call(rbind,lapply(names(child),function(v){
    k<-child[[v]];data.frame(x0=xs[as.integer(v)],x1=xs[as.integer(v)],r0=min(ys[k]),r1=max(ys[k]))
  }))
  list(lines=rbind(horizontal,vertical),x=xs,y=ys,child=child,
    tips=data.frame(isolate_id=tr$tip.label,x=xs[seq_len(nt)],row=ys[seq_len(nt)]),
    order=tr$tip.label[order],depth=max(xs[seq_len(nt)]))
}
track_data <- function(m) {
  y<-suppressWarnings(as.integer(m$year))
  year<-as.character(cut(y,c(-Inf,1999,2004,2008,2012,2016,2020,2024,Inf),
    labels=c('Before 2000','2000-2004','2005-2008','2009-2012','2013-2016','2017-2020','2021-2024','2025+')))
  year[blank(m$year)|is.na(y)]<-'Other'
  a<-m$amr_count
  amr<-as.character(cut(a,c(-Inf,0,3,6,10,13,Inf),labels=c('0','1-3','4-6','7-10','11-13','14+')))
  amr[is.na(a)]<-'Other'
  src<-ifelse(m$source=='Chatinkha nursery','Neonatal-unit sampling',m$source)
  st<-as_other(m$sequence_type);st[st!='Other']<-paste0('ST',st[st!='Other'])
  lin<-as_other(m$lineage);lin[lin!='Other']<-paste0('L',lin[lin!='Other'])
  d<-data.frame(ST=st,Lineage=lin,Country=as_other(m$country),Source=as_other(src),Year=year,
    `AMR subclass count`=amr,stringsAsFactors=FALSE,check.names=FALSE)
  allcols<-list(ST=c(ST216='#B89844',ST335='#32867B',ST5028='#8B6AB0'),
    Lineage=c(L1='#B89B44',L2='#5E977A',L3='#B46E88',L4='#5C80B2',L5='#A27D55'),
    Country=colours(d$Country),Source=colours(d$Source),
    Year=setNames(c('#EEF5DD','#D5E8CF','#ACD7C5','#79C1BB','#4FA5B2','#33869F','#286787','#214C6A'),c('Before 2000','2000-2004','2005-2008','2009-2012','2013-2016','2017-2020','2021-2024','2025+')),
    `AMR subclass count`=setNames(c('#F0EDF6','#D7CDE8','#B6A7D3','#907DB8','#6F559A','#4D3073'),c('0','1-3','4-6','7-10','11-13','14+')))
  for(k in names(d)) {
    vals<-unique(d[[k]]);missing<-setdiff(vals,c(names(allcols[[k]]),'Other'))
    if(length(missing))allcols[[k]]<-c(allcols[[k]],colours(missing))
    allcols[[k]]<-c(allcols[[k]][names(allcols[[k]])!='Other'],Other='#C2CAD0')
    allcols[[k]]<-allcols[[k]][names(allcols[[k]])%in%vals]
  }
  list(data=d,cols=allcols)
}
legend_height <- function(tk,w) {

  rows<-split(names(tk$cols),ceiling(seq_along(tk$cols)/3))
  sum(vapply(rows,function(keys)max(vapply(tk$cols[keys],length,integer(1)))*3.1+8,numeric(1)))+4
}
draw_legends <- function(tk,w,y) {
  titles<-c(ST='Sequence type',Lineage='hierBAPS lineage',Country='Country',
    Source='Source as recorded',Year='Collection year',
    `AMR subclass count`='AMR subclass count')
  xs<-9+(0:2)*(w-18)/3
  rows<-split(names(tk$cols),ceiling(seq_along(tk$cols)/3))
  for(kk in rows) {
    rowheight<-max(vapply(tk$cols[kk],length,integer(1)))*3.1+8
    for(j in seq_along(kk)) {
      k<-kk[j];col<-tk$cols[[k]];gtxt(titles[k],xs[j],y,size=7.4,bold=TRUE)
      for(i in seq_along(col)) {
        yy<-y-5-i*3.1;grect(xs[j]+1,yy,2.0,2.0,unname(col[i]))
        gtxt(names(col)[i],xs[j]+3.2,yy,size=6.4)
      }
    };y<-y-rowheight
  }
}
draw_side_legends <- function(tk,x,y,w,heading_size=6.6,item_size=5.7,pitch=3.0,gap=3.2) {
  titles<-c(ST='Sequence type',Lineage='hierBAPS lineage',Country='Country',
    Source='Source',Year='Collection year',`AMR subclass count`='AMR subclass count',
    `XDR status`='XDR status',`CARB-R status`='Carbapenem resistance',
    `MAC-R status`='Macrolide resistance',`CARB-R / MAC-R`='CARB-R / MAC-R')
  preferred<-list(
    ST=c('ST216','ST335','ST5028','Other'),
    Lineage=c('L1','L2','L3','L4','L5','Other'),
    Country=c('Malawi','Mozambique','South Africa','Uganda','United Kingdom','Nigeria',
      'Brazil','United States','Mexico','Taiwan','Other'),
    Source=c('Malawi water sample','QECH clinical isolate','Neonatal-unit sampling',
      'Eastern Cape Outbreak','Other'),
    `XDR status`=c('non-XDR','XDR','Other'),
    `CARB-R status`=c('non-CARB-R','CARB-R','Other'),
    `MAC-R status`=c('non-MAC-R','MAC-R','Other'),
    `CARB-R / MAC-R`=c('non-CARB-R / MAC-R','CARB-R / MAC-R','Other'))
  for(k in names(tk$cols)) {
    col<-tk$cols[[k]]
    if(k%in%names(preferred)) {
      wanted<-preferred[[k]];extra<-setdiff(names(col),wanted)
      col<-col[c(intersect(wanted,names(col)),extra)]
      if('Other'%in%names(col))col<-col[c(setdiff(names(col),'Other'),'Other')]
    }
    title<-if(k%in%names(titles))unname(titles[[k]]) else k
    gtxt(title,x,y,size=heading_size,bold=TRUE);y<-y-4.1
    nc<-if(length(col)>7)2L else 1L;nr<-ceiling(length(col)/nc);colw<-w/nc
    for(i in seq_along(col)) {
      cc<-ceiling(i/nr);rr<-i-(cc-1)*nr;xx<-x+(cc-1)*colw;yy<-y-(rr-1)*pitch
      if(k=='ST')grid::grid.points(xx+1.1,yy,pch=21,size=grid::unit(2.1,'mm'),
        default.units='native',gp=grid::gpar(fill=unname(col[i]),col='#44525A',lwd=.25)) else
        grect(xx+1.1,yy,2.1,2.1,unname(col[i]))
      gtxt(names(col)[i],xx+3.5,yy,size=item_size)
    }
    y<-y-nr*pitch-gap
  }
  invisible(y)
}
prepare_tree_view <- function(tr,m,binary=NULL,drop_meta=character()) {
  co<-tree_coordinates(tr);stopif(co$depth<=0,'Cannot plot a tree with zero total depth')
  m<-m[match(co$order,m$isolate_id),,drop=FALSE];stopif(anyNA(m$isolate_id),'Missing tree metadata')
  tk<-track_data(m)
  stopif(any(!drop_meta%in%names(tk$data)),'Requested metadata track to drop is absent')
  keep<-setdiff(names(tk$data),drop_meta)
  tk$data<-tk$data[keep];tk$cols<-tk$cols[keep]
  list(tr=tr,m=m,co=co,tk=tk,binary=binary)
}
rect_tree_draw <- function(view,s,labelled=TRUE,rows=NULL,page=NULL,genes=NULL) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    m<-view$m;co<-view$co;tk<-view$tk;n<-nrow(m)
    if(is.null(rows))rows<-seq_len(n)
    lo<-min(rows);hi<-max(rows);shown<-m[rows,,drop=FALSE]
    labelw<-if(labelled)max(78,min(110,w*.30)) else 0
    genenames<-genes%||%character();nmeta<-ncol(tk$data)
    right<-w-10;genew<-if(length(genenames))min(5.0,(w-115-labelw)/(length(genenames)+nmeta)) else 0
    metaw<-if(length(genenames))6.0 else 8.0
    gene_left<-right-length(genenames)*genew
    meta_left<-gene_left-nmeta*metaw-(if(length(genenames))4 else 0)
    tree_left<-10;tree_right<-meta_left-labelw-6
    stopif(tree_right-tree_left<28,'Too many tracks: reduce genes_per_block or increase canvas width')
    label_left<-tree_right+3; row_top<-h-17
    leg_h<-legend_height(tk,w);bottom<-leg_h+25
    row_denominator<-if(!is.null(page))min(n,s$cfg$rows_per_page) else length(rows)
    rowh<-(row_top-bottom)/row_denominator
    yfun<-function(r)row_top-(r-lo+.5)*rowh
    xfun<-function(x)tree_left+x/co$depth*(tree_right-tree_left)
    page_locator(w,h,page)
    if(labelled)gtxt('Isolate ID | read-run accession(s)',label_left,row_top+4,size=7,bold=TRUE)
    for(j in seq_len(nmeta))gtxt(names(tk$data)[j],meta_left+(j-.5)*metaw,row_top+3,
      size=6.5,rot=60,just='left',bold=TRUE)
    if(length(genenames)) for(j in seq_along(genenames))
      gtxt(genenames[j],gene_left+(j-.5)*genew,row_top+3,size=6.2,rot=60,just='left')

    ed<-co$lines;keep<-pmax(ed$r0,ed$r1)>=lo-.5 & pmin(ed$r0,ed$r1)<=hi+.5
    ed<-ed[keep,,drop=FALSE]
    r0<-pmax(lo-.5,pmin(hi+.5,ed$r0));r1<-pmax(lo-.5,pmin(hi+.5,ed$r1))
    if(labelled)for(r in rows[seq_along(rows)%%2==0])grect((label_left+right)/2,yfun(r),right-label_left,rowh,fill='#F5F8FA')
    gline(xfun(ed$x0),yfun(r0),xfun(ed$x1),yfun(r1),col='#52616C',lwd=if(labelled).5 else .32)
    ti<-co$tips[match(shown$isolate_id,co$tips$isolate_id),,drop=FALSE]
    if(labelled)gline(xfun(ti$x),yfun(rows),tree_right+1,yfun(rows),col='#D3DCE2',lwd=.3)
    grid::grid.points(xfun(ti$x),yfun(rows),pch=16,size=grid::unit(if(labelled).9 else .45,'mm'),
      default.units='native',gp=grid::gpar(col=unname(tk$cols$Country[tk$data$Country[rows]])))
    if(labelled)gtxt(wrap_label(shown$plot_label,labelw-4,7),label_left,yfun(rows),size=7)
    for(j in seq_len(nmeta))grect(meta_left+(j-.5)*metaw,yfun(rows),metaw*.90,rowh*.92,
      fill=unname(tk$cols[[j]][tk$data[[j]][rows]]))
    if(length(genenames)) {
      b<-view$binary
      for(j in seq_along(genenames)) {
        v<-b$present[match(paste(shown$isolate_id,genenames[j],sep='\t'),paste(b$isolate_id,b$gene,sep='\t'))]
        fill<-ifelse(is.na(v),'#BFC8CF',ifelse(v==1,'#214A64','#FFFFFF'))
        grect(gene_left+(j-.5)*genew,yfun(rows),genew*.88,rowh*.92,fill,col='#DEE5E9',lwd=.12)
      }
    }
    bar<-pretty(c(0,co$depth/4),n=3)[2];if(!is.finite(bar)||bar<=0||bar>co$depth)bar<-co$depth/5
    gline(tree_left,bottom-5,xfun(bar),bottom-5,col=INK,lwd=1)
    gtxt(paste(format(signif(bar,2)),s$cfg$branch_units),tree_left,bottom-9,size=6.7)
    if(length(genenames)) {
      gtxt('Detected',meta_left,bottom-5,size=6.7);grect(meta_left-2,bottom-5,2,2,'#214A64')
      gtxt('Not detected',meta_left+23,bottom-5,size=6.7);grect(meta_left+21,bottom-5,2,2,'white',GRID)
      gtxt('Not available',meta_left+56,bottom-5,size=6.7);grect(meta_left+54,bottom-5,2,2,'#BFC8CF')
    }
    draw_legends(tk,w,bottom-18)
    footer<-if(length(genenames))
      'Determinant and replicon tiles show recorded calls, not clinical XDR. Grey means not available, never an inferred absence.' else
      'Metadata tiles show recorded values; grey denotes Other and is not interpreted as absence.'
    gtxt(footer,9,5,size=6,col=MUTED)
  }
}
labelled_tree_page_draw <- function(view,s,rows,page,page_size) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    co<-view$co;tk<-view$tk;m<-view$m;n<-nrow(m);lo<-min(rows);hi<-max(rows)
    tree_left<-12;tree_right<-106;label_left<-112;right<-w-10
    top<-h-17;bottom<-21;rowh<-(top-bottom)/page_size
    yfun<-function(r)top-(r-lo+.5)*rowh
    xfun<-function(x)tree_left+x/co$depth*(tree_right-tree_left)
    page_locator(w,h,page)
    gtxt('Isolate ID | read-run accession(s)',label_left,top+5,size=7,bold=TRUE)
    for(r in rows[seq_along(rows)%%2==0])grect((label_left+right)/2,yfun(r),right-label_left,rowh,fill='#F5F8FA')
    ed<-co$lines;keep<-pmax(ed$r0,ed$r1)>=lo-.5&pmin(ed$r0,ed$r1)<=hi+.5
    ed<-ed[keep,,drop=FALSE]
    r0<-pmax(lo-.5,pmin(hi+.5,ed$r0));r1<-pmax(lo-.5,pmin(hi+.5,ed$r1))
    gline(xfun(ed$x0),yfun(r0),xfun(ed$x1),yfun(r1),col='#45545D',lwd=.5)
    shown<-m[rows,,drop=FALSE]
    ti<-co$tips[match(shown$isolate_id,co$tips$isolate_id),,drop=FALSE]
    gline(xfun(ti$x),yfun(rows),tree_right+1,yfun(rows),col='#D3DCE2',lwd=.3)
    country<-view$tk$data$Country[rows];cp<-view$tk$cols$Country
    grid::grid.points(xfun(ti$x),yfun(rows),pch=18,size=grid::unit(.9,'mm'),
      default.units='native',gp=grid::gpar(col=unname(cp[country])))
    gtxt(wrap_label(shown$plot_label,right-label_left-3,7),label_left,yfun(rows),size=7)
    bar<-pretty(c(0,co$depth/4),n=3)[2];if(!is.finite(bar)||bar<=0||bar>co$depth)bar<-co$depth/5
    gline(tree_left,12,xfun(bar),12,col=INK,lwd=1)
    gtxt(format(signif(bar,2)),tree_left,8,size=6.5)
  }
}
manuscript_tree_page_draw <- function(view,s,rows,page,page_size) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    co<-view$co;m<-view$m;lo<-min(rows);hi<-max(rows)
    tree_left<-8;tree_right<-68;label_left<-73;right<-w-8
    top<-h-15;bottom<-20;rowh<-(top-bottom)/page_size
    yfun<-function(r)top-(r-lo+.5)*rowh
    xfun<-function(x)tree_left+x/co$depth*(tree_right-tree_left)
    page_locator(w,h,page)
    gtxt('Isolate ID | read-run accession(s)',label_left,top+5,size=7,bold=TRUE)
    for(r in rows[seq_along(rows)%%2==0])grect((label_left+right)/2,yfun(r),right-label_left,rowh,
      fill='#F5F8FA')
    ed<-co$lines;keep<-pmax(ed$r0,ed$r1)>=lo-.5&pmin(ed$r0,ed$r1)<=hi+.5
    ed<-ed[keep,,drop=FALSE]
    r0<-pmax(lo-.5,pmin(hi+.5,ed$r0));r1<-pmax(lo-.5,pmin(hi+.5,ed$r1))
    gline(xfun(ed$x0),yfun(r0),xfun(ed$x1),yfun(r1),col='#39474F',lwd=.55)
    shown<-m[rows,,drop=FALSE]
    ti<-co$tips[match(shown$isolate_id,co$tips$isolate_id),,drop=FALSE]
    gline(xfun(ti$x),yfun(rows),tree_right+1,yfun(rows),col='#D3DCE2',lwd=.3)
    country<-view$tk$data$Country[rows];cp<-view$tk$cols$Country
    grid::grid.points(xfun(ti$x),yfun(rows),pch=18,size=grid::unit(1.0,'mm'),
      default.units='native',gp=grid::gpar(col=unname(cp[country])))
    gtxt(wrap_label(shown$plot_label,right-label_left-3,7),label_left,yfun(rows),size=7)
    bar<-pretty(c(0,co$depth/4),n=3)[2];if(!is.finite(bar)||bar<=0||bar>co$depth)bar<-co$depth/5
    gline(tree_left,11,xfun(bar),11,col=INK,lwd=1)
    gtxt(format(signif(bar,2)),tree_left,7,size=6.5)
  }
}
write_isolate_tracking_pages <- function(s,view,prefix) {
  nr<-nrow(view$m);page_size<-min(40L,as.integer(s$cfg$rows_per_page))
  starts<-seq.int(1,nr,by=page_size);membership<-list()
  for(pg in seq_along(starts)) {
    rows<-seq.int(starts[pg],min(nr,starts[pg]+page_size-1L))
    locator<-sprintf('Page %d/%d | Isolates %d-%d/%d',pg,length(starts),min(rows),max(rows),nr)
    file<-file.path(s$out,'supplementary',paste0(prefix,'_PAGE_',sprintf('%02d',pg),'.tiff'))
    save_tiff_drawing(s,file,labelled_tree_page_draw(view,s,rows,locator,page_size),297,210)
    manuscript_file<-file.path(s$out,'supplementary',paste0(prefix,
      '_MANUSCRIPT_PAGE_',sprintf('%02d',pg),'.tiff'))
    save_tiff_drawing(s,manuscript_file,
      manuscript_tree_page_draw(view,s,rows,locator,page_size),180,240)
    membership[[pg]]<-data.frame(page=pg,display_order=rows,isolate_id=view$co$order[rows])
  }
  membership<-do.call(rbind,membership)
  stopif(!identical(membership$isolate_id,view$co$order),'Compact isolate-tracking pagination failed')
  export(s,membership,paste0(prefix,'_page_membership'))
  caption(s,prefix,paste(
    'Two isolate-tracking series show every tree tip once in fixed phylogenetic order: A4-landscape PAGE files and portrait-manuscript MANUSCRIPT_PAGE files.',
    'Labels give the analysis isolate ID and read-run accession; page membership and the plotted Newick are supplied as source data.',
    'The 180 x 240 mm MANUSCRIPT_PAGE TIFFs are intended for direct placement at full text width on consecutive portrait supplementary pages; they carry no embedded title, so one overarching external figure caption can serve the complete series.'))
  writeLines(c(
    paste0('SUPPLEMENT INSERTION GUIDE: ',prefix),
    'Preferred for a portrait supplementary document: insert MANUSCRIPT_PAGE TIFFs in numerical order, one image per page, at 180 mm wide (or the maximum text width without enlargement).',
    'Keep all pages together under one supplementary figure number; put the full caption before page 1 and use a short continuation line in the document header/footer if required by the journal.',
    'Alternative for a landscape section: use the PAGE TIFFs in numerical order at full landscape-page width.',
    'Do not crop the page locator, scale bar, tree or isolate/run labels.',
    paste0('Machine-readable coverage: source_data/',prefix,'_page_membership.csv')),
    file.path(s$out,'source_data',paste0(prefix,'_INSERTION_GUIDE.txt')))
  invisible(membership)
}
write_constituent_tree_pages <- function(s,view,specs,prefix) {
  need('ape');manifest<-list();membership<-list();full_order<-view$co$order
  if(!length(specs))return(invisible(list(manifest=data.frame(),membership=list())))
  for(i in seq_along(specs)) {
    spec<-specs[[i]];ids<-intersect(full_order,unique(spec$ids));slug<-safe_slug(spec$label)
    if(length(ids)<2L) {
      manifest[[i]]<-data.frame(group=spec$label,rule=spec$rule,n=length(ids),pages=0,
        state='SKIPPED: fewer than two verified isolates',stringsAsFactors=FALSE)
      next
    }
    tr<-ape::keep.tip(view$tr,ids)
    m<-view$m[match(tr$tip.label,view$m$isolate_id),,drop=FALSE]
    subview<-prepare_tree_view(tr,m,drop_meta=character())
    inherited<-full_order[full_order%in%ids]
    if(!identical(subview$co$order,inherited))note(s,paste(
      'Pruning changed the traversal order for',spec$label,
      '; both subgroup order and inherited master positions are exported explicitly.'))
    page_size<-min(40L,as.integer(s$cfg$rows_per_page));starts<-seq.int(1,length(ids),by=page_size)
    for(pg in seq_along(starts)) {
      rows<-seq.int(starts[pg],min(length(ids),starts[pg]+page_size-1L))
      locator<-sprintf('%s | Page %d/%d | Isolates %d-%d/%d',spec$label,pg,length(starts),
        min(rows),max(rows),length(ids))
      file<-file.path(s$out,'supplementary',paste0(prefix,'_',slug,
        '_MANUSCRIPT_PAGE_',sprintf('%02d',pg),'.tiff'))
      save_tiff_drawing(s,file,manuscript_tree_page_draw(subview,s,rows,locator,page_size),180,240)
      membership[[length(membership)+1L]]<-data.frame(group=spec$label,rule=spec$rule,
        page=pg,subgroup_order=rows,inherited_tree_order=match(subview$co$order[rows],full_order),
        isolate_id=subview$co$order[rows],stringsAsFactors=FALSE)
    }
    ape::write.tree(tr,file=file.path(s$out,'source_data',paste0(prefix,'_',slug,'_pruned_from_master.nwk')))
    export(s,subview$m,paste0(prefix,'_',slug,'_tip_metadata'))
    manifest[[i]]<-data.frame(group=spec$label,rule=spec$rule,n=length(ids),pages=length(starts),
      state='BUILT',stringsAsFactors=FALSE)
  }
  manifest<-do.call(rbind,manifest);export(s,manifest,paste0(prefix,'_manifest'))
  if(length(membership))export(s,do.call(rbind,membership),paste0(prefix,'_membership'))
  caption(s,prefix,paste(
    'Portrait 180 x 240 mm isolate-labelled constituent trees are pruned views of the same supplied master topology, never separately inferred trees.',
    'Each biological group is split into at most 40 tips per page, labels include isolate ID and read-run accession, and inherited master-tree positions are supplied in the membership table.',
    'Pages contain no embedded figure title; keep each group together under one external supplementary caption.'))
  writeLines(c(paste('SUPPLEMENT CONSTITUENT-TREE INSERTION GUIDE:',prefix),
    'Use the 180 x 240 mm MANUSCRIPT_PAGE TIFFs at full text width, in group then page-number order.',
    'Each group is a pruned view of the same master topology. Do not describe these as independently inferred trees.',
    paste0('Group definitions and page coverage: source_data/',prefix,'_manifest.csv and ',
      prefix,'_membership.csv.')),
    file.path(s$out,'source_data',paste0(prefix,'_INSERTION_GUIDE.txt')))
  invisible(list(manifest=manifest,membership=membership))
}
write_tree_outputs <- function(s,tr,m,prefix,title,binary=NULL,main=FALSE,drop_meta=character()) {
  need('ape');v<-prepare_tree_view(tr,m,binary,drop_meta)
  allgenes<-if(is.null(binary))character() else unique(binary$gene)
  blocksize<-as.integer(s$cfg$genes_per_block);blocks<-if(length(allgenes))split(allgenes,ceiling(seq_along(allgenes)/blocksize)) else list(character())
  export(s,v$m,paste0(prefix,'_tip_metadata'))
  export(s,data.frame(display_order=seq_along(v$co$order),isolate_id=v$co$order),paste0(prefix,'_tip_order'))
  export(s,data.frame(isolate_id=v$m$isolate_id,v$tk$data,check.names=FALSE),paste0(prefix,'_display_annotations'))
  export(s,v$co$lines,paste0(prefix,'_fixed_branch_coordinates'))
  ape::write.tree(tr,file=file.path(s$out,'source_data',paste0(prefix,'_plotted.nwk')))
  if(!is.null(binary))export(s,binary,paste0(prefix,'_determinant_states'))
  pages<-list();nr<-nrow(m)
  for(b in seq_along(blocks)) {
    gs<-blocks[[b]];suffix<-if(length(blocks)>1)paste0('_genes',sprintf('%02d',b)) else ''
    width<-if(length(gs))420 else 297
    leg<-legend_height(v$tk,width)
    labelw<-max(78,min(110,width*.30))
    estimated_lines<-max(1,ceiling(max(nchar(v$m$plot_label))/(floor((labelw-4)/(7*.18)))))
    pitch<-max(3.6,2.8*estimated_lines)
    fullh<-18+nr*pitch+leg+25
    starts<-seq.int(1,nr,by=as.integer(s$cfg$rows_per_page))
    ph<-max(210,18+min(nr,s$cfg$rows_per_page)*pitch+leg+25)
    if(isTRUE(s$cfg$vector_archive)) {

      archive<-file.path(s$out,'vector_archive','supplementary',paste0(prefix,suffix,'_ARCHIVAL_ALL_IDS.pdf'))
      save_pdf_drawing(archive,rect_tree_draw(v,s,TRUE,genes=gs),width,fullh)
      grDevices::pdf(file.path(s$out,'vector_archive','supplementary',paste0(prefix,suffix,'_READABLE_PAGES.pdf')),
        width=width/25.4,height=ph/25.4,useDingbats=FALSE)
      tryCatch({for(pg in seq_along(starts)) {
        rows<-seq.int(starts[pg],min(nr,starts[pg]+s$cfg$rows_per_page-1))
        page<-sprintf('Page %d/%d | Tips %d-%d/%d',pg,length(starts),min(rows),max(rows),nr)
        grid::grid.newpage();rect_tree_draw(v,s,TRUE,rows=rows,page=page,genes=gs)(width,ph)
      }},finally=grDevices::dev.off())
    }
    for(pg in seq_along(starts)) {
      rows<-seq.int(starts[pg],min(nr,starts[pg]+s$cfg$rows_per_page-1))
      page<-sprintf('Page %d/%d | Tips %d-%d/%d',pg,length(starts),min(rows),max(rows),nr)
      page_file<-file.path(s$out,'supplementary',paste0(prefix,suffix,'_READABLE_PAGE_',sprintf('%02d',pg),'.tiff'))
      save_tiff_drawing(s,page_file,rect_tree_draw(v,s,TRUE,rows=rows,page=page,genes=gs),width,ph)
      pages[[length(pages)+1]]<-data.frame(gene_block=b,page=pg,display_order=rows,isolate_id=v$co$order[rows])
    }
    if(main) {
      draw<-rect_tree_draw(v,s,FALSE,genes=gs)
      save_drawing(s,paste0('Figure_4_ST335',suffix),draw,w=max(180,95+length(gs)*4.5),h=18+nr*.47+leg+25)
    }
  }
  membership<-do.call(rbind,pages);export(s,membership,paste0(prefix,'_page_membership'))
  for(b in seq_along(blocks))stopif(!identical(membership$isolate_id[membership$gene_block==b],v$co$order),
    'Pagination failed: each isolate must appear exactly once per gene block, in full-tree order')
  caption(s,prefix,paste(title,
    'All labelled outputs use one fixed topology, one x scale and one tip order. Pages are windows of that tree, not separately reconstructed subtrees.',
    'Numbered READABLE_PAGE TIFFs contain every tip in fixed full-tree order.',
    if(isTRUE(s$cfg$vector_archive))'Optional combined READABLE_PAGES and one-page ARCHIVAL_ALL_IDS PDFs are retained in vector_archive.' else '',
    'Local labels are analysis isolate IDs followed by run accessions; public genomes already named by a run accession are shown once. Unresolved mappings are explicitly marked.',
    'Columns describe recorded metadata. Replicon detection does not establish plasmid co-location; determinant calls are not phenotypes.'))
  invisible(v)
}

circular_draw <- function(view,s) {
  function(w,h) {
    stopif(w<=h,'CC25 main figure must use a landscape canvas')
    stopif(w/h<1.22||w/h>1.40,'CC25 canvas aspect ratio must remain between 1.22 and 1.40')
    co<-view$co;tk<-view$tk;n<-nrow(view$m)
    ring_keys<-setdiff(names(tk$data),'ST')
    stopif(length(ring_keys)!=5L,
      'CC25 main figure requires exactly five rings: lineage, country, source, year and AMR count')
    legendw<-56;treew<-w-legendw
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2,
      widths=grid::unit(c(treew,legendw),'mm'))))
    grid::pushViewport(grid::viewport(layout.pos.col=1,xscale=c(0,treew),
      yscale=c(0,h),clip='on'))
    scale_band<-15;outer<-min((treew-8)/2,(h-scale_band-6)/2)
    ringw<-2.1;R<-outer-length(ring_keys)*ringw-3
    stopif(R<32,'Circular tree canvas is too small')
    cx<-treew/2;cy<-scale_band+3+outer

    gap<-.17;theta_start<-pi/2-gap/2
    theta<-function(row)theta_start-(row-.5)/n*(2*pi-gap)
    rad<-function(x)x/co$depth*R
    for(i in seq_len(nrow(view$tr$edge))) {
      a<-view$tr$edge[i,1];b<-view$tr$edge[i,2];ang<-theta(co$y[b])
      gline(cx+rad(co$x[a])*cos(ang),cy+rad(co$x[a])*sin(ang),
        cx+rad(co$x[b])*cos(ang),cy+rad(co$x[b])*sin(ang),col='#56656F',lwd=.38)
    }
    for(k in names(co$child)) {
      children<-co$child[[k]];a<-theta(range(co$y[children]));r<-rad(co$x[as.integer(k)])
      th<-seq(a[1],a[2],length.out=max(3,ceiling(abs(diff(a))*50)))
      grid::grid.lines(cx+r*cos(th),cy+r*sin(th),default.units='native',gp=grid::gpar(col='#56656F',lwd=.38))
    }
    ti<-co$tips[match(view$m$isolate_id,co$tips$isolate_id),,drop=FALSE]
    stcol<-unname(tk$cols$ST[tk$data$ST])
    grid::grid.points(cx+rad(ti$x)*cos(theta(seq_len(n))),cy+rad(ti$x)*sin(theta(seq_len(n))),
      pch=21,size=grid::unit(.90,'mm'),default.units='native',
      gp=grid::gpar(fill=stcol,col='#44525A',lwd=.25))

    for(j in seq_along(ring_keys)) {
      rin<-R+2+(j-1)*ringw;rout<-rin+ringw*.9
      for(i in seq_len(n)) {
        th<-seq(theta(i-.46),theta(i+.46),length.out=4)
        xx<-cx+c(rin*cos(th),rout*cos(rev(th)));yy<-cy+c(rin*sin(th),rout*sin(rev(th)))
        grid::grid.polygon(xx,yy,default.units='native',
          gp=grid::gpar(fill=unname(tk$cols[[ring_keys[j]]][tk$data[[ring_keys[j]]][i]]),col=NA))
      }
    }
    short<-c(Lineage='Lineage',Country='Country',Source='Source',Year='Year',
      `AMR subclass count`='AMR count')
    for(j in seq_along(ring_keys)) {
      rmid<-R+2+(j-.5)*ringw
      gtxt(unname(short[[ring_keys[j]]]),cx,cy+rmid,size=5.1,just='centre',bold=TRUE)
    }
    bar<-pretty(c(0,co$depth/4),n=3)[2];if(!is.finite(bar)||bar<=0||bar>co$depth)bar<-co$depth/5
    yscale<-8;xscale<-max(8,cx-rad(bar)/2)
    gline(xscale,yscale,xscale+rad(bar),yscale,col=INK,lwd=1)
    gtxt(format(signif(bar,2)),xscale+rad(bar)/2,yscale-4,size=6.2,just='centre')
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.col=2,xscale=c(0,legendw),
      yscale=c(0,h),clip='on'))
    legend_bottom<-draw_side_legends(tk,3,h-7,legendw-6)
    stopif(legend_bottom<4,'CC25 legend overflow: increase the right legend panel height')
    grid::popViewport();grid::popViewport()
  }
}
module_cc25 <- function(s) {
  ob<-load_tree_set(s)
  v<-write_tree_outputs(s,ob$tree,ob$meta,'Supplement_CC25','CC25: unique isolate and read-run labels')
  write_isolate_tracking_pages(s,v,'Supplement_CC25_isolate_tracking')
  lineage<-as_other(v$m$lineage);specs<-list()
  for(g in sort(setdiff(unique(lineage),'Other'))) {
    display<-if(grepl('^L',g,ignore.case=TRUE))toupper(g) else paste0('L',g)
    specs[[length(specs)+1L]]<-list(label=paste0('Lineage_',display),
      rule=paste('CC25 hierBAPS lineage',display),ids=v$m$isolate_id[lineage==g])
  }
  write_constituent_tree_pages(s,v,specs,'Supplement_CC25_constituent_trees')
  draw<-circular_draw(v,s);save_drawing(s,'Figure_3A_CC25',draw,210,165)
  caption(s,'Figure_3A_CC25',paste(
    'CC25 maximum-likelihood tree, displayed using the supplied branch lengths and a midpoint root when configured.',
    'Tip-circle colour denotes sequence type. From inside to outside the five annotation rings show hierBAPS lineage, country, source, collection year and recorded AMR-subclass count.',
    'Short headers sit directly in the narrow top seam, while compact, vertically spaced keys occupy a physically separate right-hand legend panel.',
    'Lineage is encoded as a ring rather than unverified monophyletic-clade shading. Every tip is retained.',
    'Unavailable or uncategorised metadata are labelled Other. Supplement_CC25 contains full annotation pages, and Supplement_CC25_isolate_tracking contains compact A4-landscape pages with isolate/read-run labels in identical tip order.',
    'Branch lengths are those of the midpoint-rooted maximum-likelihood tree inferred from the recombination-filtered SNP alignment and are shown in model branch-length units; they are not rescaled into SNP counts.'))

  if(has(s,'master')||(has(s,'genomes')&&has(s,'metadata'))) {
    tryCatch({master<-load_master(s);outside<-master[!master$isolate_id%in%ob$tree$tip.label,,drop=FALSE]
      export(s,outside,'genomes_outside_CC25_tree')
    },error=function(e)note(s,paste('Full-genome register not completed:',conditionMessage(e))))
  }
}
load_binary <- function(s,m) {
  if(has(s,'amr')) {
    a<-read_input(s,'amr');require_columns(a,c('isolate_id','gene','present'),'AMR calls')
    a$isolate_id<-clean_id(a$isolate_id);a$gene<-trimws(a$gene)
    unique_ids(paste(a$isolate_id,a$gene,sep='\t'),'AMR isolate/gene pairs')
    stopif(any(blank(a$gene)),'AMR gene names are missing')
    a$present<-numeric_checked(a$present,'AMR present')
    stopif(any(!is.na(a$present)&!a$present%in%c(0,1)),'AMR calls must be 0, 1 or blank/not available')
    export(s,a[!a$isolate_id%in%m$isolate_id,,drop=FALSE],'AMR_rows_outside_ST335_plot')
    grid<-expand.grid(isolate_id=m$isolate_id,gene=unique(a$gene),stringsAsFactors=FALSE)
    grid$present<-a$present[match(paste(grid$isolate_id,grid$gene,sep='\t'),paste(a$isolate_id,a$gene,sep='\t'))]
    a<-grid
  } else {
    master<-load_master(s)
    possible<-grep('^(bla|qnr|aac|aph|aad|arr|mph|msr|mdf|gyrA|gyrB|parC|parE|floR|cat|cml|erm|fos|sul|dfr|tet|qac)',names(master),value=TRUE)
    cols<-possible[vapply(master[possible],function(x)all(blank(x)|x%in%c(0,1)),logical(1))]
    stopif(!length(cols),'No explicit 0/1 AMR gene columns in master. Supply publication_inputs/amr_determinants.csv.')
    stopif(!all(m$isolate_id%in%master$isolate_id),'ST335 tip missing from master AMR table')
    a<-do.call(rbind,lapply(cols,function(g)data.frame(isolate_id=m$isolate_id,gene=g,
      present=numeric_checked(master[[g]][match(m$isolate_id,master$isolate_id)],g))))
    note(s,'AMR calls taken from explicit binary gene columns of the master table, not XDR categories.')
  }
  stopif(any(a$gene%in%c('IncC','IncHI2')),'Keep IncC/IncHI2 in the replicon table, not duplicate gene columns')
  if(has(s,'replicons')) {
    r<-identify(read_input(s,'replicons'));require_columns(r,c('IncC','IncHI2'),'Replicon table')
    r<-do.call(rbind,lapply(c('IncC','IncHI2'),function(k)data.frame(isolate_id=m$isolate_id,gene=k,
      present=numeric_checked(r[[k]][match(m$isolate_id,r$isolate_id)],k))))
    stopif(any(!is.na(r$present)&!r$present%in%c(0,1)),'Replicon calls must be 0, 1 or blank/not available')
  } else {
    master<-tryCatch(load_master(s),error=function(e)NULL)
    available<-if(is.null(master))character() else intersect(c('IncC','IncHI2'),names(master))
    r<-do.call(rbind,lapply(c('IncC','IncHI2'),function(k){
      value<-rep(NA_real_,nrow(m))
      if(k%in%available) {
        ii<-match(m$isolate_id,master$isolate_id)
        value<-numeric_checked(master[[k]][ii],k)
        stopif(any(!is.na(value)&!value%in%c(0,1)),paste(k,'master-table calls must be 0, 1 or blank/not available'))
      }
      data.frame(isolate_id=m$isolate_id,gene=k,present=value)
    }))
    missing_replicons<-setdiff(c('IncC','IncHI2'),available)
    if(!length(available)) {
      if(isTRUE(s$cfg$strict))stop('Figure 4 requires verified IncC/IncHI2 calls in strict mode')
      note(s,'Replicon input missing: IncC and IncHI2 are shown as not available, not absent.')
    } else {
      if(length(missing_replicons)) {
        if(isTRUE(s$cfg$strict))stop(paste('Missing verified master-table replicon call:',
          paste(missing_replicons,collapse=', ')))
        note(s,paste('Replicon calls unavailable for',paste(missing_replicons,collapse=', '),
          '; those cells are shown as not available, not absent.'))
      }
      note(s,paste('Replicon calls were taken from explicit master-table columns:',
        paste(available,collapse=', ')))
    }
  }
  a<-rbind(a,r)

  if(anyNA(a$present))note(s,paste(sum(is.na(a$present)),'AMR/replicon cells are not available in the plotted subset.'))
  a
}
get_st335 <- function(s) {
  if(!is.null(s$st335))return(s$st335)
  ob<-load_tree_set(s);ids<-ob$meta$isolate_id[!blank(ob$meta$sequence_type)&ob$meta$sequence_type=='335']
  stopif(length(ids)<2,'Fewer than two ST335 tips have verified MLST')
  dedicated<-s$cfg$st335_tree %||% ''
  if(nzchar(dedicated)) {
    p<-if(grepl('^(/|[A-Za-z]:)',dedicated))dedicated else file.path(s$root,dedicated)
    tr<-ape::read.tree(p);s$used<-unique(c(s$used,p));stopif(!inherits(tr,'phylo'),'ST335 input must be one tree')
    tr$tip.label<-clean_id(tr$tip.label);unique_ids(tr$tip.label,'ST335 tree tips')
    stopif(!setequal(tr$tip.label,ids),'Dedicated ST335 tree tip set differs from ST335 membership in CC25')
    stopif(is.null(tr$edge.length)||any(!is.finite(tr$edge.length))||any(tr$edge.length<0),'Invalid ST335 branch lengths')
    if(isTRUE(s$cfg$midpoint_root))tr<-phangorn::midpoint(tr)
    tr<-ape::ladderize(tr)
    origin<-'Supplied ST335-specific Newick.'
  } else {
    tr<-ape::keep.tip(ob$tree,ids)
    origin<-'ST335-only view pruned from the supplied CC25 tree; not a newly inferred or independently fitted ST335 phylogeny.'
    note(s,origin)
  }
  m<-ob$meta[match(tr$tip.label,ob$meta$isolate_id),,drop=FALSE]
  if(nrow(m)!=224)note(s,paste('ST335 has',nrow(m),'tips; manuscript reports 224.'))
  s$st335<-list(tree=tr,meta=m,origin=origin);s$st335
}
status_from_values <- function(x,positive,negative,pos_label,neg_label) {
  z<-tolower(trimws(as.character(x)));out<-rep('Other',length(z))
  out[!blank(z)&z%in%tolower(positive)]<-pos_label
  out[!blank(z)&z%in%tolower(negative)]<-neg_label
  out
}
st335_status_tracks <- function(s,m,binary=NULL) {
  unique_ids(m$isolate_id,'ST335 status metadata IDs')
  sources<-list(`CC25 metadata`=m)
  if(!is.null(s$master)||has(s,'master')||(has(s,'genomes')&&has(s,'metadata'))) {
    master<-tryCatch(load_master(s),error=function(e){
      note(s,paste('Master-table status fields could not be loaded; per-isolate unavailable values are retained:',
        conditionMessage(e)));NULL
    })
    if(!is.null(master))sources[['master table']]<-master
  }
  status_missing<-function(x) {
    z<-tolower(trimws(as.character(x)))
    blank(x)|z%in%c('unknown','unavailable','not available','not known','other','na','n/a')
  }
  status_token<-function(x) {
    z<-tolower(trimws(as.character(x)));out<-z
    out[z%in%c('1','true','yes','positive','resistant','xdr','carb-r','mac-r',
      'carb-r / mac-r','carb-r/mac-r','carb+r mac+r','both')]<-'positive'
    out[z%in%c('0','false','no','negative','susceptible','non-xdr','non xdr','not xdr',
      'non-carb-r','non-mac-r','non-carb-r / mac-r','non-carb-r/mac-r','neither')]<-'negative'
    out
  }
  recorded_field<-function(candidates) {
    value<-rep(NA_character_,nrow(m));field_by_row<-source_by_row<-rep('',nrow(m))
    for(nm in candidates)for(src_name in names(sources)) {
      src<-sources[[src_name]];if(!nm%in%names(src))next
      ii<-match(m$isolate_id,src$isolate_id);incoming<-rep(NA_character_,nrow(m))
      ok<-!is.na(ii);incoming[ok]<-as.character(src[[nm]][ii[ok]])
      incoming[status_missing(incoming)]<-NA_character_
      overlap<-!status_missing(value)&!status_missing(incoming)
      conflict<-overlap&status_token(value)!=status_token(incoming)
      stopif(any(conflict),paste('Conflicting recorded status values for',
        paste(head(m$isolate_id[conflict],8),collapse=', '),
        '- reconcile',paste(unique(c(field_by_row[conflict],nm)),collapse=' / ')))
      fill<-status_missing(value)&!status_missing(incoming)
      value[fill]<-incoming[fill];field_by_row[fill]<-nm;source_by_row[fill]<-src_name
    }
    value[status_missing(value)]<-NA_character_;used<-!status_missing(value)
    list(value=value,
      field=if(any(used))paste(unique(field_by_row[used]),collapse='; ') else '',
      source=if(any(used))paste(unique(source_by_row[used]),collapse='; ') else '',
      field_by_row=field_by_row,source_by_row=source_by_row)
  }
  xfield<-recorded_field(c('XDR (Traditional)','XDR status','XDR_status','XDR.Traditional','XDR'))
  xraw<-xfield$value
  xdr<-status_from_values(xraw,c('1','true','yes','xdr','resistant'),
    c('0','false','no','non-xdr','non xdr','not xdr','susceptible'),'XDR','non-XDR')
  if(!nzchar(xfield$field))note(s,'No recorded XDR status column was found; the XDR track is shown as Other.')
  carbfield<-recorded_field(c('CARB-R status','CARB-R','CARB_R','carbapenem resistance',
    'Carbapenem resistance','carbapenem_status'))
  macfield<-recorded_field(c('MAC-R status','MAC-R','MAC_R','macrolide resistance',
    'Macrolide resistance','macrolide_status'))
  combinedfield<-recorded_field(c('CARB-R / MAC-R','CARB-R/MAC-R','CARB_MAC_R','CARB.MAC.R',
    'CARB-R / MAC-R status','CARB_MAC_status'))
  carb_status<-status_from_values(carbfield$value,c('1','true','yes','positive','carb-r','resistant'),
    c('0','false','no','negative','non-carb-r','susceptible'),'CARB-R','non-CARB-R')
  mac_status<-status_from_values(macfield$value,c('1','true','yes','positive','mac-r','resistant'),
    c('0','false','no','negative','non-mac-r','susceptible'),'MAC-R','non-MAC-R')
  craw<-combinedfield$value

  combined<-status_from_values(craw,
    c('1','true','yes','positive','carb-r / mac-r','carb-r/mac-r','carb+r mac+r','both'),
    c('0','false','no','negative','non-carb-r / mac-r','non-carb-r/mac-r','neither'),
    'CARB-R / MAC-R','non-CARB-R / MAC-R')
  n<-nrow(m);carb_detected<-mac_detected<-additional_mac_detected<-rep('',n)
  carb_complete<-mac_complete<-additional_mac_complete<-rep(FALSE,n)
  marker_source<-''
  combined_carb_status<-combined_mac_status<-rep('Other',n)
  marker_component<-function(markers,pos_label,neg_label) {
    state<-rep('Other',n);detected<-rep('',n);complete<-rep(FALSE,n)
    if(is.null(binary)||!length(markers))return(list(state=state,detected=detected,complete=complete))
    keys<-paste(binary$isolate_id,binary$gene,sep='\t')
    for(i in seq_len(n)) {
      v<-binary$present[match(paste(m$isolate_id[i],markers,sep='\t'),keys)]
      hit<-markers[!is.na(v)&v==1];detected[i]<-paste(hit,collapse=';')
      complete[i]<-length(v)>0&&all(!is.na(v))
      state[i]<-if(length(hit))pos_label else if(complete[i])neg_label else 'Other'
    }
    list(state=state,detected=detected,complete=complete)
  }
  if(!is.null(binary)) {
    gs<-unique(binary$gene)
    gene_key<-toupper(gsub('[^A-Za-z0-9]','',gs))
    carb<-gs[grepl('^(BLANDM|BLAKPC|BLAVIM|BLAIMP|BLAOXA48)',gene_key)]
    mac<-gs[grepl('^(MPH|MSR|ERM)',gene_key)]

    additional_mac<-gs[gene_key%in%c('MPHE','MSRE')]
    cs<-marker_component(carb,'CARB-R','non-CARB-R')
    ms<-marker_component(mac,'MAC-R','non-MAC-R')
    ams<-marker_component(additional_mac,'MAC-R','non-MAC-R')
    carb_detected<-cs$detected;mac_detected<-ms$detected
    additional_mac_detected<-ams$detected
    carb_complete<-cs$complete;mac_complete<-ms$complete;additional_mac_complete<-ams$complete
    fill_carb<-blank(carbfield$value);carb_status[fill_carb]<-cs$state[fill_carb]
    fill_mac<-blank(macfield$value);mac_status[fill_mac]<-ms$state[fill_mac]

    combined_carb_status<-cs$state;combined_mac_status<-ams$state
    if(length(carb)&&length(additional_mac))marker_source<-paste(
      'For isolates lacking an explicit combined field, CARB-R / MAC-R requires both a recorded carbapenemase-family call and an additional mphE/msrE call; mphA alone is not sufficient.',
      'Complete non-positive component calls are non-CARB-R / MAC-R and incomplete component calls are Other.')
  }
  derived<-ifelse(combined_carb_status=='CARB-R'&combined_mac_status=='MAC-R','CARB-R / MAC-R',
    ifelse(combined_carb_status!='Other'&combined_mac_status!='Other','non-CARB-R / MAC-R','Other'))
  use_derived<-blank(craw);combined[use_derived]<-derived[use_derived]
  if(!nzchar(combinedfield$field)&&!nzchar(marker_source)) {
    note(s,'No explicit combined field and no complete carbapenemase plus mphE/msrE evidence were available; CARB-R / MAC-R is shown as Other where it cannot be classified.')
  }
  evidence<-data.frame(isolate_id=m$isolate_id,xdr_source_value=xraw,
    xdr_source_field=xfield$field_by_row,xdr_source_table=xfield$source_by_row,
    carb_source_value=carbfield$value,carb_source_field=carbfield$field_by_row,
    carb_source_table=carbfield$source_by_row,mac_source_value=macfield$value,
    mac_source_field=macfield$field_by_row,mac_source_table=macfield$source_by_row,
    combined_source_value=craw,combined_source_field=combinedfield$field_by_row,
    combined_source_table=combinedfield$source_by_row,
    carbapenem_markers_detected=carb_detected,
    macrolide_markers_detected=mac_detected,
    additional_mphE_msrE_detected=additional_mac_detected,
    carbapenem_calls_complete=carb_complete,macrolide_calls_complete=mac_complete,
    additional_mphE_msrE_calls_complete=additional_mac_complete,
    `CARB-R status`=carb_status,`MAC-R status`=mac_status,
    combined_carbapenem_marker_component=combined_carb_status,
    combined_additional_mphE_msrE_component=combined_mac_status,
    combined_derivation_rule=ifelse(!blank(craw),
      paste('Recorded field',combinedfield$field_by_row,'from',combinedfield$source_by_row),marker_source),
    check.names=FALSE,stringsAsFactors=FALSE)
  out<-data.frame(isolate_id=m$isolate_id,`XDR status`=xdr,
    `CARB-R / MAC-R`=combined,check.names=FALSE,stringsAsFactors=FALSE)
  export(s,cbind(out,evidence[setdiff(names(evidence),'isolate_id')]),'Figure_4A_ST335_status_tracks')
  counts<-do.call(rbind,lapply(setdiff(names(out),'isolate_id'),function(k)
    data.frame(track=k,state=names(table(out[[k]])),n=as.integer(table(out[[k]])))))
  export(s,counts,'Figure_4A_ST335_status_counts')
  if(n==224L) {
    nx<-sum(xdr=='XDR');nc<-sum(combined=='CARB-R / MAC-R')
    if(nx!=199L)note(s,paste('Recorded XDR count is',nx,'of 224; the manuscript reports 199. Inspect the exported status audit.'))
    if(nc!=5L)note(s,paste('Combined CARB-R / MAC-R count is',nc,'of 224; the manuscript reports 5. Inspect the exported status audit.'))
  }
  attr(out,'component_status')<-data.frame(isolate_id=m$isolate_id,
    `CARB-R status`=carb_status,`MAC-R status`=mac_status,check.names=FALSE,
    stringsAsFactors=FALSE)
  attr(out,'methods')<-c(
    if(nzchar(xfield$field))paste('XDR status uses recorded field',xfield$field,'from',xfield$source,'.') else
      'XDR status was unavailable and is displayed as Other.',
    if(nzchar(combinedfield$field))paste('CARB-R / MAC-R uses recorded field',combinedfield$field,
      'from',combinedfield$source,'where populated;',marker_source) else if(nzchar(marker_source))marker_source else
      'CARB-R / MAC-R status was unavailable and is displayed as Other.')
  out
}
st335_landscape_draw <- function(view,s,genes=character(),status=NULL) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    co<-view$co;m<-view$m;n<-nrow(m);top<-h-10;bottom<-30
    stopif(is.null(status),'ST335 main view requires auditable XDR and CARB-R / MAC-R tracks')
    unique_ids(status$isolate_id,'ST335 status IDs')
    status<-status[match(m$isolate_id,status$isolate_id),,drop=FALSE]
    stopif(anyNA(status$isolate_id),'ST335 status rows do not cover every plotted isolate')
    status_names<-c('XDR status','CARB-R / MAC-R')
    stopif(any(!status_names%in%names(status)),'ST335 status table is missing a required display column')
    tree_left<-8;tree_right<-74;year_left<-tree_right+7;year_w<-5.5
    legend_w<-58;legend_left<-w-legend_w+3
    status_w<-5.8;status_right<-legend_left-7
    status_left<-status_right-length(status_names)*status_w
    gene_left<-year_left+year_w+3;gene_right<-status_left-4
    gene_w<-if(length(genes))(gene_right-gene_left)/length(genes) else 0
    stopif(length(genes)>0&&gene_w<2.25,'Too many determinant columns for the fixed paper-ready ST335 canvas')
    rowh<-(top-bottom)/n
    yfun<-function(r)top-(r-.5)*rowh
    xfun<-function(x)tree_left+x/co$depth*(tree_right-tree_left)
    ed<-co$lines
    gline(xfun(ed$x0),yfun(ed$r0),xfun(ed$x1),yfun(ed$r1),col='#242C31',lwd=.68)
    ti<-co$tips[match(m$isolate_id,co$tips$isolate_id),,drop=FALSE]
    country<-as_other(m$country)
    country_base<-c(Malawi='#E31A1C',Mozambique='#377EB8',`South Africa`='#4DAF4A',
      Uganda='#984EA3',`United Kingdom`='#F28E2B',Other='#B8BEC2')
    extra_country<-setdiff(unique(country),names(country_base))
    country_cols<-c(country_base,colours(extra_country))
    country_cols<-country_cols[!duplicated(names(country_cols))]
    country_cols<-country_cols[names(country_cols)%in%unique(country)]
    if('Other'%in%names(country_cols))country_cols<-country_cols[c(setdiff(names(country_cols),'Other'),'Other')]
    grid::grid.points(xfun(ti$x),yfun(seq_len(n)),pch=18,
      size=grid::unit(min(1.10,max(.95,rowh*1.65)),'mm'),default.units='native',
      gp=grid::gpar(col=unname(country_cols[country]),fill=unname(country_cols[country]),lwd=.25))
    years<-as_other(m$year);known<-sort(unique(years[years!='Other']))
    year_cols<-if(length(known))setNames(grDevices::colorRampPalette(
      c('#F3EFA6','#9BD5C2','#3A9CA6','#173D7A'))(length(known)),known) else character()
    year_cols<-c(year_cols,Other='#C3C8CB')
    grect(year_left+year_w/2,yfun(seq_len(n)),year_w*.92,rowh*.94,
      fill=unname(year_cols[years]),col=NA)
    gtxt('Year',year_left+year_w/2,bottom-3,size=6.5,rot=55,just='right',bold=TRUE)
    if(length(genes)) {
      b<-view$binary
      keys<-paste(b$isolate_id,b$gene,sep='\t')
      for(j in seq_along(genes)) {
        v<-b$present[match(paste(m$isolate_id,genes[j],sep='\t'),keys)]
        fill<-ifelse(is.na(v),'#B8BEC2',ifelse(v==1,'#1696A7','#F4F6F7'))
        grect(gene_left+(j-.5)*gene_w,yfun(seq_len(n)),gene_w*.88,rowh*.94,
          fill=fill,col='#E3E7E9',lwd=.08)
        gtxt(genes[j],gene_left+(j-.5)*gene_w,bottom-3,size=6.2,rot=58,just='right')
      }
    }
    for(j in seq_along(status_names)) {
      vals<-status[[status_names[j]]]
      positive<-vals%in%c('XDR','CARB-R / MAC-R')
      negative<-vals%in%c('non-XDR','non-CARB-R / MAC-R')
      fill<-ifelse(positive,'#F05A61',ifelse(negative,'#D7D9DB','#FFFFFF'))
      grect(status_left+(j-.5)*status_w,yfun(seq_len(n)),status_w*.86,rowh*.94,
        fill=fill,col='#D6DADD',lwd=.08)
      gtxt(status_names[j],status_left+(j-.5)*status_w,bottom-3,size=6.2,rot=58,just='right')
    }
    bar<-pretty(c(0,co$depth/4),n=3)[2];if(!is.finite(bar)||bar<=0||bar>co$depth)bar<-co$depth/5
    gline(tree_left,12,xfun(bar),12,col=INK,lwd=1)
    gtxt(format(signif(bar,2)),tree_left,8,size=6.4)
    ly<-top-2
    gtxt('Country',legend_left,ly,size=7.8,bold=TRUE)
    for(i in seq_along(country_cols)) {
      yy<-ly-4-i*3.35
      grid::grid.points(legend_left+1,yy,pch=18,size=grid::unit(1.05,'mm'),default.units='native',
        gp=grid::gpar(col=unname(country_cols[i])))
      gtxt(names(country_cols)[i],legend_left+4.2,yy,size=6.8)
    }
    ly<-ly-7-length(country_cols)*3.35
    gtxt('Collection year',legend_left,ly,size=7.8,bold=TRUE)
    year_show<-year_cols[names(year_cols)%in%unique(years)]
    for(i in seq_along(year_show)) {
      yy<-ly-4-i*3.05;grect(legend_left+1,yy,2.2,2.2,unname(year_show[i]))
      gtxt(names(year_show)[i],legend_left+4.2,yy,size=6.6)
    }
    if(length(genes)) {
      ly<-ly-7-length(year_show)*3.05
      gtxt('AMR determinant',legend_left,ly,size=7.8,bold=TRUE)
      fills<-c(Present='#1696A7',Absent='#F4F6F7',`Not available`='#B8BEC2')
      for(i in seq_along(fills)) {
        yy<-ly-4-i*3.25;grect(legend_left+1,yy,2.2,2.2,unname(fills[i]),'#D2D6D8',.2)
        gtxt(names(fills)[i],legend_left+4.2,yy,size=6.6)
      }
      ly<-ly-8-length(fills)*3.25
    }
    status_keys<-list(`XDR status`=c(`non-XDR`='#D7D9DB',XDR='#F05A61',Other='#FFFFFF'),
      `CARB-R / MAC-R`=c(`non-CARB-R / MAC-R`='#D7D9DB',`CARB-R / MAC-R`='#F05A61',Other='#FFFFFF'))
    for(k in status_names) {
      gtxt(k,legend_left,ly,size=7.8,bold=TRUE);cols<-status_keys[[k]]
      for(i in seq_along(cols)) {
        yy<-ly-4-i*3.25;grect(legend_left+1,yy,2.2,2.2,unname(cols[i]),'#D2D6D8',.2)
        gtxt(names(cols)[i],legend_left+4.2,yy,size=6.6)
      }
      ly<-ly-8-length(cols)*3.25
    }
    stopif(ly<3,'ST335 legend overflow: reduce displayed exact-year levels or determinant columns')
  }
}
curate_st335_genes <- function(s,a,requested='') {
  if(is.null(a))return(character())
  asked<-trimws(strsplit(requested%||%'',',',fixed=TRUE)[[1]]);asked<-asked[nzchar(asked)]
  genes<-unique(a$gene)
  if(length(asked)) {
    stopif(any(!asked%in%genes),'A requested main_genes name is absent from the call table; check exact spelling')
    stopif(length(unique(asked))>30L,'The fixed paper-ready ST335 panel supports at most 30 requested determinant columns')
    return(unique(asked))
  }
  recorded<-tapply(a$present,a$gene,function(x)any(!is.na(x)))
  available<-names(recorded)[recorded]
  key<-setNames(toupper(gsub('[^A-Za-z0-9]','',available)),available)
  targets<-list(
    `blaCTX-M-15`=c('^BLACTXM15'),`blaOXA-48`=c('^BLAOXA48'),
    `blaNDM-1`=c('^BLANDM1'),`blaTEM-1`=c('^BLATEM1'),
    `aac(6)-Ib`=c('^AAC6IB$'),`aac(6)-Ib-cr5`=c('^AAC6IBCR5'),
    `aac(6)-IIe`=c('^AAC6IIE'),
    `aac(3)-IId`=c('^AAC3IID$','^AAC3II$'),`aac(3)-IIe`=c('^AAC3IIE'),
    aadA1=c('^AADA1'),aadA2=c('^AADA2'),`aph(3)-Ib`=c('^APH3IB'),
    `aph(6)-Id`=c('^APH6ID'),arr=c('^ARR'),qnrB1=c('^QNRB1'),
    sul1=c('^SUL1'),sul2=c('^SUL2'),`mph(A)`=c('^MPHA$'),`mph(E)`=c('^MPHE$'),
    `msr(E)`=c('^MSRE$'),`tet(A)`=c('^TETA$'),cmlA5=c('^CMLA5'),catA1=c('^CATA1'),
    floR=c('^FLOR'),dfrA12=c('^DFRA12'),dfrA14=c('^DFRA14'),dfrA23=c('^DFRA23'),
    qacEdelta1=c('^QACEDELTA1','^QACE1'),IncC=c('^INCC$'),IncHI2=c('^INCHI2$'))
  selected<-character();mapping<-data.frame(target=names(targets),selected_gene=NA_character_)
  for(i in seq_along(targets)) {
    hit<-character()
    for(pattern in targets[[i]]) {
      hit<-setdiff(names(key)[grepl(pattern,key)],selected)
      if(length(hit))break
    }
    if(length(hit)) {
      chosen<-sort(hit)[1];selected<-c(selected,chosen);mapping$selected_gene[i]<-chosen
    }
  }
  selected<-unique(selected)
  if(!length(selected)) {
    positive<-tapply(a$present,a$gene,function(x)any(!is.na(x)&x==1))
    selected<-head(names(positive)[positive],25L)
    note(s,'No curated manuscript determinant names matched; the main ST335 panel uses the first 25 recorded-positive determinants and the inventory documents this fallback.')
  }
  attr(selected,'target_mapping')<-mapping
  selected
}
module_st335_tree <- function(s) {
  ob<-get_st335(s)
  binary_source<-has(s,'amr')||!is.null(s$master)||has(s,'master')||
    (has(s,'genomes')&&has(s,'metadata'))
  a<-if(binary_source)load_binary(s,ob$meta) else {
    note(s,'No determinant source was supplied; determinant columns are unavailable rather than inferred.')
    NULL
  }
  v<-write_tree_outputs(s,ob$tree,ob$meta,'Supplement_ST335_labelled',
    if(is.null(a))'ST335: unique isolate and read-run labels with metadata' else
      'ST335: unique isolate and read-run labels with metadata and recorded determinant calls',
    binary=a,drop_meta='ST')
  status<-st335_status_tracks(s,v$m,a)
  write_isolate_tracking_pages(s,v,'Supplement_ST335_isolate_tracking')
  country_group<-as_other(v$m$country);lineage_group<-as_other(v$m$lineage)
  source_group<-as_other(ifelse(v$m$source=='Chatinkha nursery','Neonatal-unit sampling',v$m$source))
  constituent_specs<-list(
    list(label='Malawi',rule='ST335 isolates recorded in Malawi',
      ids=v$m$isolate_id[country_group=='Malawi']),
    list(label='South_Africa',rule='ST335 isolates recorded in South Africa',
      ids=v$m$isolate_id[country_group=='South Africa']))
  for(g in c('4','5','L4','L5'))if(g%in%lineage_group)constituent_specs[[length(constituent_specs)+1L]]<-
    list(label=paste0('South_Africa_lineage_',safe_slug(g)),
      rule=paste('South African ST335 in hierBAPS lineage',g),
      ids=v$m$isolate_id[country_group=='South Africa'&lineage_group==g])
  residual<-country_group!='Other'&!country_group%in%c('Malawi','South Africa')
  if(any(residual))constituent_specs[[length(constituent_specs)+1L]]<-
    list(label='Other_countries',rule='ST335 isolates with a recorded country outside Malawi and South Africa',
      ids=v$m$isolate_id[residual])
  eastern<-grepl('eastern cape',tolower(source_group))
  if(any(eastern))constituent_specs[[length(constituent_specs)+1L]]<-
    list(label='Eastern_Cape_outbreak',rule='Source recorded as Eastern Cape outbreak',
      ids=v$m$isolate_id[eastern])
  for(g in intersect(c('Malawi water sample','QECH clinical isolate','Neonatal-unit sampling'),
    unique(source_group)))constituent_specs[[length(constituent_specs)+1L]]<-
    list(label=paste0('Malawi_source_',safe_slug(g)),rule=paste('Malawi ST335 source recorded as',g),
      ids=v$m$isolate_id[country_group=='Malawi'&source_group==g])
  if(sum(status[['XDR status']]=='non-XDR')>=2L)constituent_specs[[length(constituent_specs)+1L]]<-
    list(label='non_XDR',rule='Recorded genotypic XDR status is non-XDR',
      ids=status$isolate_id[status[['XDR status']]=='non-XDR'])
  if(sum(status[['CARB-R / MAC-R']]=='CARB-R / MAC-R')>=2L)
    constituent_specs[[length(constituent_specs)+1L]]<-
      list(label='CARB_R_MAC_R',rule='Combined CARB-R / MAC-R status is positive',
        ids=status$isolate_id[status[['CARB-R / MAC-R']]=='CARB-R / MAC-R'])
  write_constituent_tree_pages(s,v,constituent_specs,'Supplement_ST335_constituent_trees')
  selected<-curate_st335_genes(s,a,s$cfg$main_genes)
  target_mapping<-attr(selected,'target_mapping')
  if(!is.null(target_mapping))export(s,target_mapping,'Figure_4A_ST335_curated_target_mapping')
  inventory<-if(is.null(a))data.frame(gene=character(),has_recorded_call=logical(),in_main=logical()) else {
    has<-tapply(a$present,a$gene,function(x)any(!is.na(x)&x==1))
    recorded<-tapply(a$present,a$gene,function(x)any(!is.na(x)))
    data.frame(gene=names(has),has_recorded_call=as.logical(recorded[names(has)]),
      has_positive_call=as.logical(has),in_main=names(has)%in%selected,
      main_order=match(names(has),selected))
  }
  export(s,inventory,'Figure_4A_ST335_gene_inventory')
  save_drawing(s,'Figure_4A_ST335_tree_heatmap',st335_landscape_draw(v,s,selected,status),250,160)
  caption(s,'Supplement_ST335_labelled',paste(ob$origin,
    'The invariant sequence-type column is omitted because every displayed isolate is ST335.',
    'Numbered 600-dpi READABLE_PAGE TIFFs contain every displayed isolate in fixed tree order; compact A4 isolate-tracking pages are also supplied for direct placement in supplementary material.',
    if(isTRUE(s$cfg$vector_archive))'Optional combined and full-height PDFs are retained in vector_archive.' else '',
    'Local isolate labels include the matching run accession; public genomes already named by a run accession are shown once.',
    'Metadata tiles show recorded values; grey denotes Other and is not interpreted as absence.',
    if(is.null(a))'No determinant table was available for this render.' else
      'Determinant tiles show recorded calls: dark teal present, near-white absent and grey not available.'))
  caption(s,'Figure_4A_ST335_tree_heatmap',paste(ob$origin,
    'Landscape ST335 phylogeny with enlarged country-coloured tip diamonds, exact collection-year strip, aligned recorded AMR-determinant columns and separated right-hand legends.',
    if(is.null(a))'Determinant columns could not be drawn because no valid determinant input was available; the tree, country and year remain shown.' else
      'The main panel uses a fixed manuscript-led determinant order, or the exact main_genes selection when supplied; the complete call table remains in the labelled supplement.',
    paste(attr(status,'methods'),collapse=' '),
    'Unavailable or uncategorised metadata are labelled Other. Branch lengths remain in model branch-length units and are not converted into SNP counts.'))
}
module_snp_matrices <- function(s) {
  cc<-load_tree_set(s);cc_order<-tree_coordinates(cc$tree)$order
  cc_meta<-cc$meta[match(cc_order,cc$meta$isolate_id),,drop=FALSE]
  cc_out<-write_snp_outputs(s,cc_order,'CC25',require_exact=TRUE,meta=cc_meta,
    audit_extras=TRUE)
  st<-get_st335(s);st_order<-tree_coordinates(st$tree)$order
  st_meta<-st$meta[match(st_order,st$meta$isolate_id),,drop=FALSE]
  a<-if(has(s,'amr'))load_binary(s,st_meta) else tryCatch(load_binary(s,st_meta),error=function(e){
    note(s,paste('SNP metadata status components could not use determinant calls:',conditionMessage(e)));NULL
  })
  status<-st335_status_tracks(s,st_meta,a)
  st_out<-write_snp_outputs(s,st_order,'ST335',require_exact=FALSE,meta=st_meta,
    status=status,audit_extras=FALSE)
  manifest<-list(
    data.frame(panel='CC25_full',rule='All isolates in the plotted CC25 tree',n=nrow(cc_out$D),
      labelled=FALSE,state='BUILT',filename='Supplement_CC25_pairwise_SNP_heatmap.tiff'),
    data.frame(panel='ST335_full',rule='All isolates in the plotted ST335 tree',n=nrow(st_out$D),
      labelled=FALSE,state='BUILT',filename='Supplement_ST335_pairwise_SNP_heatmap.tiff'))
  membership<-list();add_result<-function(z) {
    manifest[[length(manifest)+1L]]<<-z$manifest
    if(nrow(z$membership))membership[[length(membership)+1L]]<<-z$membership
    invisible(NULL)
  }

  st_value<-as_other(cc_meta$sequence_type)
  add_result(write_snp_panel(s,cc_out$D,cc_meta,cc_meta$isolate_id[st_value=='216'],
    'CC25_ST216','Sequence type recorded as ST216','CC25',common_max=max(cc_out$D)))
  lineage<-as_other(cc_meta$lineage)
  for(g in sort(setdiff(unique(lineage),'Other')))if(sum(lineage==g)>=2L)
    add_result(write_snp_panel(s,cc_out$D,cc_meta,cc_meta$isolate_id[lineage==g],
      paste0('CC25_lineage_',safe_slug(g)),paste('hierBAPS lineage',g),'CC25',
      common_max=max(cc_out$D)))
  country<-as_other(st_meta$country)
  for(g in intersect(c('Malawi','South Africa'),unique(country)))
    add_result(write_snp_panel(s,st_out$D,st_meta,st_meta$isolate_id[country==g],
      paste0('ST335_country_',safe_slug(g)),paste('ST335 country recorded as',g),'ST335',status,
      max(st_out$D)))
  st_lineage<-as_other(st_meta$lineage)
  for(g in intersect(c('4','5','L4','L5'),unique(st_lineage))) {
    ids<-st_meta$isolate_id[country=='South Africa'&st_lineage==g]
    if(length(ids)>=2L)add_result(write_snp_panel(s,st_out$D,st_meta,ids,
      paste0('ST335_South_Africa_lineage_',safe_slug(g)),
      paste('South African ST335 in hierBAPS lineage',g),'ST335',status,max(st_out$D)))
  }
  source<-as_other(ifelse(st_meta$source=='Chatinkha nursery','Neonatal-unit sampling',st_meta$source))
  for(g in intersect(c('Malawi water sample','QECH clinical isolate','Neonatal-unit sampling'),
    unique(source))) {
    ids<-st_meta$isolate_id[country=='Malawi'&source==g]
    if(length(ids)>=2L)add_result(write_snp_panel(s,st_out$D,st_meta,ids,
      paste0('ST335_Malawi_source_',safe_slug(g)),paste('Malawi ST335 source recorded as',g),
      'ST335',status,max(st_out$D)))
  }
  eastern<-grepl('eastern cape',tolower(source))
  if(sum(eastern)>=2L)add_result(write_snp_panel(s,st_out$D,st_meta,st_meta$isolate_id[eastern],
    'ST335_Eastern_Cape_outbreak','ST335 source identifies the Eastern Cape outbreak',
    'ST335',status,max(st_out$D)))
  for(spec in list(c('XDR status','non-XDR'),c('CARB-R / MAC-R','CARB-R / MAC-R'))) {
    k<-spec[1];g<-spec[2]
    if(sum(status[[k]]==g)>=2L)add_result(write_snp_panel(s,st_out$D,st_meta,
      status$isolate_id[status[[k]]==g],paste0('ST335_',safe_slug(k),'_',safe_slug(g)),
      paste(k,'recorded or reproducibly derived as',g),'ST335',status,max(st_out$D)))
  }
  manifest<-do.call(rbind,manifest);export(s,manifest,'SNP_panel_manifest')
  if(length(membership))export(s,do.call(rbind,membership),'SNP_subgroup_membership')
  st_group<-as_other(cc_meta$sequence_type);st_group[st_group!='Other']<-paste0('ST',st_group[st_group!='Other'])
  lineage_group<-as_other(cc_meta$lineage);lineage_group[lineage_group!='Other']<-paste0('L',lineage_group[lineage_group!='Other'])
  comparisons<-list(
    pairwise_group_comparison(cc_out$D,st_group,
      'CC25 sequence type'),
    pairwise_group_comparison(cc_out$D,lineage_group,
      'CC25 hierBAPS lineage'),
    pairwise_group_comparison(st_out$D,country,'ST335 country'),
    pairwise_group_comparison(st_out$D,source,'ST335 source'),
    pairwise_group_comparison(st_out$D,status[['XDR status']],'ST335 XDR status'),
    pairwise_group_comparison(st_out$D,status[['CARB-R / MAC-R']],'ST335 CARB-R / MAC-R status'))
  comparison_long<-do.call(rbind,lapply(comparisons,'[[','long'))
  comparison_summary<-do.call(rbind,lapply(comparisons,'[[','summary'))
  if(nrow(comparison_long))export(s,comparison_long,'SNP_pairwise_distance_comparisons_long')
  if(nrow(comparison_summary))export(s,comparison_summary,'SNP_pairwise_distance_comparisons_summary')
  caption(s,'Supplement_CC25_pairwise_SNP_heatmap',paste(
    'Pairwise SNP-distance matrix for the exact CC25 tree-tip set, ordered identically to the plotted phylogeny.',
    'Mirrored top and left strips show sequence type, lineage, country, source, collection-year band and AMR-subclass count for every isolate.',
    'The numeric matrix is supplied in wide TSV and CSV formats and as a unique-pair long table. The colour transform is square-root scaled but its key is labelled in raw SNP counts.'))
  caption(s,'Supplement_ST335_pairwise_SNP_heatmap',paste(
    'Exact ST335 subset of the validated CC25 pairwise SNP-distance matrix, ordered identically to the plotted ST335 view.',
    'Mirrored metadata strips show lineage, country, source, year, XDR status and combined CARB-R / MAC-R status for every isolate.',
    'Distances are directly subset from the supplied SNP matrix; they are not inferred from phylogenetic branch lengths. The colour transform is square-root scaled but its key is labelled in raw SNP counts.'))
  writeLines(c(
    'NUMERIC SNP MATRICES',
    'CC25_pairwise_SNP_matrix_tree_order.tsv and .csv',
    'ST335_pairwise_SNP_matrix_tree_order.tsv and .csv',
    'The corresponding *_pairwise_SNP_long.csv, *_pairwise_SNP_summary.csv and *_SNP_isolate_order_metadata.csv files are in this same source_data directory.',
    'SNP_panel_manifest.csv lists the full and biologically stratified panels.',
    'SNP_subgroup_membership.csv records every inclusion and inherited matrix position.',
    'SNP_pairwise_distance_comparisons_long.csv and _summary.csv provide descriptive within/between-group comparisons without non-independent pairwise p-values.',
    'Rendered 600-dpi TIFF heatmaps with mirrored metadata strips and right-hand legends are in the sibling supplementary directory.'),
    file.path(s$out,'source_data','README_SNP_MATRICES.txt'))
}
module_st335 <- function(s) {
  ob<-get_st335(s);tr<-ob$tree;m<-ob$meta;origin<-ob$origin
  a<-load_binary(s,m)

  write_tree_outputs(s,tr,m,'Supplement_ST335','ST335: unique isolate labels and determinant calls',a,FALSE,drop_meta='ST')
  selected<-trimws(strsplit(s$cfg$main_genes%||%'',',',fixed=TRUE)[[1]])
  selected<-selected[nzchar(selected)]
  if(!length(selected))selected<-unique(a$gene)
  stopif(any(!selected%in%a$gene),'A requested main_genes name is absent from the call table; check exact spelling')
  selected<-unique(c(setdiff(selected,c('IncC','IncHI2')),'IncC','IncHI2'))
  export(s,data.frame(gene=unique(a$gene),in_main=unique(a$gene)%in%selected),'Figure_4_gene_inventory')
  v<-prepare_tree_view(tr,m,a,drop_meta='ST')
  blocks<-split(selected,ceiling(seq_along(selected)/as.integer(s$cfg$genes_per_block)))
  if(length(blocks)>1)note(s,paste('Main Figure 4 split into',length(blocks),'gene blocks; select main_genes explicitly for one compact panel. All genes remain in supplements.'))
  for(b in seq_along(blocks)) {
    gs<-blocks[[b]];width<-max(200,110+length(gs)*5.0)
    h<-18+nrow(m)*.50+legend_height(v$tk,width)+22
    name<-if(length(blocks)==1)'Figure_4B_ST335_determinants' else paste0('Figure_4B_ST335_determinants_part',sprintf('%02d',b))
    save_drawing(s,name,rect_tree_draw(v,s,FALSE,genes=gs),width,h)
    caption(s,name,paste(origin,'Individual gene calls and replicon calls are displayed separately from clinical phenotypes.',
      'Dark: detected; white: not detected in a completed analysis; grey: not available.',
      'No genotypic XDR, CARB-R or MAC-R classification is assigned. Raw values and all tip labels accompany the figure.'))
  }
}

KEN_ALL <- c('AMP10','SXT25','C30','PEF5','CIP5','CPD10','CTX5','AZM15','MEM10','IPM10',
  'FOX30','FEP30','TCG15','AK30','ATM30','CN10','TZP36')
ken_draw <- function(z,m,agents,all=FALSE) {
  function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    n<-nrow(z);left<-if(all)93 else 76;right<-w-10;top<-h-17;bottom<-19
    rowh<-(top-bottom)/n;cw<-(right-left)/length(agents)
    y<-top-(seq_len(n)-.5)*rowh
    labels<-paste(m$isolate_id,m$ena_run,sep=' | ')
    gtxt('Isolate ID | read-run accession',9,top+4,size=7.2,bold=TRUE)
    for(i in seq_len(n))if(i%%2==0)grect((9+right)/2,y[i],right-9,rowh,fill='#F4F7F9')
    gtxt(labels,9,y,size=7.8)
    ramp<-grDevices::colorRampPalette(c('#EEF4F6','#9EC7CD','#3A8F9B','#184B67'))(101)
    for(j in seq_along(agents)) {
      g<-agents[j];v<-numeric_checked(z[[g]],g);x<-left+(j-.5)*cw
      idx<-pmax(1,pmin(101,round((v-6)/30*100)+1))
      fill<-ifelse(is.na(v),'#C2CAD0',ramp[idx])
      grect(x,y,cw*.95,rowh*.92,fill=fill)
      gtxt(ifelse(is.na(v),'NA',v),x,y,size=7.8,
        just='centre',col=ifelse(is.na(v)|v<25,INK,'#FFFFFF'))
      gtxt(g,x,top+4,size=6.8,just='left',rot=55,bold=TRUE)
    }
    by<-9;bx<-left;bw<-min(44,right-left-30)
    for(i in 1:100)grect(bx+(i-.5)/100*bw,by,bw/100+.02,2.6,fill=ramp[i])
    gtxt('6',bx,by-3,size=6.3);gtxt('36 mm',bx+bw,by-3,size=6.3,just='right')
    gtxt('Recorded zone diameter',bx,by+3.8,size=6.5)
  }
}
module_ken <- function(s) {
  z<-read_bundled(s,'ken_zones.csv');m<-read_bundled(s,'ken_metadata.csv')
  require_columns(z,c('isolate_id',KEN_ALL),'Ken zones');unique_ids(z$isolate_id,'Ken zone IDs');unique_ids(m$isolate_id,'Ken metadata IDs')
  stopif(!setequal(z$isolate_id,m$isolate_id),'Ken zone and metadata ID sets differ')
  m<-m[match(z$isolate_id,m$isolate_id),,drop=FALSE]
  for(g in KEN_ALL) {
    z[[g]]<-numeric_checked(z[[g]],g)
    stopif(any(!is.na(z[[g]])&(z[[g]]<6|z[[g]]>100|z[[g]]!=round(z[[g]]))),paste('Invalid raw zone:',g))
  }

  mm<-attach_accessions(s,identify(m));stopif(any(canonical_runs(m$ena_run)!=mm$run_accession),'Ken accession mismatch')
  save_drawing(s,'Figure_2_Ken_zones',ken_draw(z,m,KEN_ALL,TRUE),297,210)
  save_drawing(s,'Supplement_Ken_all_17_antibiotics',ken_draw(z,m,KEN_ALL,TRUE),297,210,TRUE)
  export(s,merge(z,m,by='isolate_id',sort=FALSE),'Figure_2_all_zones_and_metadata')
  export(s,data.frame(display_order=seq_len(nrow(z)),isolate_id=z$isolate_id),'Figure_2_row_order')
  key<-data.frame(code=KEN_ALL,antibiotic=c('Ampicillin','Co-trimoxazole','Chloramphenicol','Pefloxacin',
    'Ciprofloxacin','Cefpodoxime','Cefotaxime','Azithromycin','Meropenem','Imipenem',
    'Cefoxitin','Cefepime','Tigecycline','Amikacin','Aztreonam','Gentamicin','Piperacillin-tazobactam'),
    disk_content_microgram=c('10','1.25/23.75','30','5','5','10','5','15','10','10',
      '30','30','15','30','30','10','30/6'))
  export(s,key,'Figure_2_drug_key')
  caption(s,'Figure_2_Ken_zones',paste('Recorded disk-diffusion zones for the 19-isolate neonatal-unit subset.',
    'Cell colour and values show zone diameter in mm, not MIC or a universal resistance category.',
    'Rows retain workbook order. Labels give isolate name and the matching run in the supplied accession workbook.',
    'No composite XDR label is plotted. All 17 tested antibiotics are retained in the main and supplementary heatmaps, and the drug key accompanies both.'))
}

load_malawi <- function(s) {
  z<-read_input(s,'malawi');require_columns(z,c('Isolate','Class'),'Malawi metadata/matrix')
  ids<-clean_id(z$Isolate);unique_ids(ids,'Malawi isolates')
  names_clean<-clean_id(names(z));ix<-which(names_clean%in%ids)
  stopif(length(ix)!=length(ids)||anyDuplicated(names_clean[ix])>0,
    'Cannot match SNP columns to all Malawi isolate IDs. Export named columns; no positional column guessing is used.')
  D<-as.matrix(data.frame(lapply(z[ix],numeric_checked,what='Malawi SNP counts'),check.names=FALSE))
  rownames(D)<-ids;colnames(D)<-names_clean[ix];D<-validate_distance(D)
  m<-z[,-ix,drop=FALSE];m$isolate_id<-ids
  m$source<-as_other(ifelse(m$Class=='Chatinkha nursery','Neonatal-unit sampling',m$Class))
  m<-attach_accessions(s,m)
  list(meta=m,D=D)
}
kruskal <- function(D) {
  ii<-which(upper.tri(D),arr.ind=TRUE)
  e<-data.frame(a=rownames(D)[ii[,1]],b=colnames(D)[ii[,2]],snps=D[ii])
  e<-e[order(e$snps,e$a,e$b),,drop=FALSE];comp<-setNames(seq_len(nrow(D)),rownames(D));keep<-logical(nrow(e));count<-0
  if(nrow(e))for(i in seq_len(nrow(e))) {
    a<-comp[e$a[i]];b<-comp[e$b[i]]
    if(a!=b) {keep[i]<-TRUE;comp[comp==b]<-a;count<-count+1;if(count==nrow(D)-1)break}
  }
  e<-e[keep,,drop=FALSE];stopif(nrow(D)>1&&nrow(e)!=nrow(D)-1,'MST is not connected');e
}
collapse_zero <- function(D) {

  ids<-sort(rownames(D));D<-D[ids,ids,drop=FALSE];group<-setNames(rep(NA_integer_,length(ids)),ids);k<-0L
  for(id in ids)if(is.na(group[id])) {
    k<-k+1L;near<-ids[D[id,]==0&is.na(group)]
    stopif(any(D[near,near,drop=FALSE]!=0),'Zero-SNP grouping is not transitive; use collapse_zero=false')
    stopif(any(vapply(near,function(z)any(D[z,]!=D[id,]),logical(1))),
      'Zero-SNP isolates have different pairwise profiles; possible variable callable masks. Use collapse_zero=false.')
    group[near]<-k
  }
  membership<-data.frame(isolate_id=ids,node_id=sprintf('G%03d',group[ids]),stringsAsFactors=FALSE)
  reps<-vapply(split(ids,group[ids]),function(v)v[1],character(1))
  gd<-D[reps,reps,drop=FALSE];rownames(gd)<-colnames(gd)<-sprintf('G%03d',seq_along(reps))
  list(D=gd,membership=membership)
}
network_plot <- function(nodes,edges,membership,palette) {
  need(c('ggplot2','ggrepel'))
  e<-edges
  for(k in c('x','y')) {
    e[[k]]<-nodes[[k]][match(e$a,nodes$node_id)]
    e[[paste0(k,'end')]]<-nodes[[k]][match(e$b,nodes$node_id)]
  }
  span<-max(diff(range(nodes$x)),diff(range(nodes$y)),1)
  radius<-sqrt(nodes$n/max(nodes$n))*span*.07
  slices<-list()
  for(i in seq_len(nrow(nodes))) {
    fractions<-table(membership$source[membership$node_id==nodes$node_id[i]])
    ends<-c(0,cumsum(as.numeric(fractions))/sum(fractions))*2*pi
    for(j in seq_along(fractions)) {
      th<-seq(ends[j],ends[j+1],length.out=max(5,ceiling((ends[j+1]-ends[j])*30)))
      slices[[length(slices)+1]]<-data.frame(
        x=nodes$x[i]+c(0,cos(th),0)*radius[i],
        y=nodes$y[i]+c(0,sin(th),0)*radius[i],
        slice=paste(i,j,sep='_'),source=names(fractions)[j])
    }
  }
  pieces<-do.call(rbind,slices)
  nodes$label<-paste0(nodes$node_id,' (',nodes$n,')')
  ggplot2::ggplot()+ggplot2::geom_segment(data=e,ggplot2::aes(x=x,y=y,xend=xend,yend=yend),
      colour='#AEBBC4',linewidth=.45)+
    ggplot2::geom_label(data=e,ggplot2::aes(x=(x+xend)/2,y=(y+yend)/2,label=snps),
      size=2.3,label.padding=grid::unit(.10,'lines'),fill='white',colour=MUTED)+
    ggplot2::geom_polygon(data=pieces,ggplot2::aes(x=x,y=y,group=slice,fill=source),colour='white',linewidth=.2)+
    ggplot2::scale_fill_manual(values=palette,name='Source')+
    ggrepel::geom_text_repel(data=nodes,ggplot2::aes(x=x,y=y,label=label),seed=20260919L,
      size=2.6,min.segment.length=0,max.overlaps=Inf,box.padding=.5,
      max.time=Inf,max.iter=20000,point.padding=.6)+
    ggplot2::coord_equal(clip='off')+ggplot2::theme_void(base_family='sans')+
    ggplot2::theme(legend.position='bottom',
      legend.title=ggplot2::element_text(size=8),legend.text=ggplot2::element_text(size=7),
      plot.margin=ggplot2::margin(7,10,7,10))
}
module_network <- function(s) {
  need(c('ggplot2','ggrepel','igraph'));ob<-load_malawi(s);D<-ob$D;m<-ob$meta
  if(isTRUE(s$cfg$collapse_zero))g<-collapse_zero(D) else {
    ids<-sort(rownames(D));g<-list(D=D[ids,ids,drop=FALSE],membership=data.frame(isolate_id=ids,node_id=ids))
  }
  membership<-merge(g$membership,m[,c('isolate_id','run_accession','plot_label','source')],by='isolate_id',sort=FALSE)
  gd<-g$D;edges<-kruskal(gd);nid<-rownames(gd)
  if(length(nid)==1)xy<-matrix(c(0,0),1) else {
    graph<-igraph::graph_from_data_frame(edges[,c('a','b')],directed=FALSE,vertices=data.frame(name=nid))
    set.seed(s$cfg$seed);xy<-igraph::layout_with_fr(graph,weights=1/(1+edges$snps),niter=2000)
  }
  nodes<-data.frame(node_id=nid,x=xy[,1],y=xy[,2],stringsAsFactors=FALSE)
  nodes$n<-vapply(nid,function(g)sum(membership$node_id==g),integer(1))
  nodes$source<-vapply(nid,function(g){v<-unique(membership$source[membership$node_id==g]);if(length(v)==1)v else 'Mixed sources'},character(1))
  nodes$label<-nodes$node_id
  pal<-colours(membership$source)
  p<-network_plot(nodes,edges,membership,pal)
  save_gg(s,p,'Figure_1A_Malawi_SNP_network',180,150)
  export(s,as.data.frame(table(membership$node_id,membership$source)),'Figure_1A_node_source_counts');
  export(s,membership,'Figure_1A_node_isolate_run_membership');export(s,nodes,'Figure_1A_node_coordinates');export(s,edges,'Figure_1A_SNP_edges')

  mem<-membership[order(membership$node_id,membership$isolate_id),,drop=FALSE]
  draw<-function(w,h) {
    native(w,h);on.exit(grid::popViewport())
    y<-h-14-(seq_len(nrow(mem))-.5)*3.7
    for(i in seq_len(nrow(mem)))if(i%%2==0)grect(w/2,y[i],w-18,3.7,'#F4F7F9')
    gtxt('Node',9,h-9,7,bold=TRUE);gtxt('Isolate ID | read-run accession',35,h-9,7,bold=TRUE);gtxt('Source',168,h-9,7,bold=TRUE)
    gtxt(mem$node_id,9,y,7);gtxt(mem$plot_label,35,y,7);gtxt(mem$source,168,y,7)
  }
  save_drawing(s,'Supplement_Malawi_network_all_IDs',draw,297,20+nrow(mem)*3.7,supp=TRUE,pdf=FALSE)
  export_snp_matrix(s,D,'Malawi')
  matrix_meta<-m[match(rownames(D),m$isolate_id),,drop=FALSE]
  stopif(anyNA(matrix_meta$isolate_id),'Malawi SNP heatmap metadata failed exact matrix-order alignment')
  source_value<-as_other(matrix_meta$source)
  source_cols<-colours(source_value)
  source_cols<-c(source_cols[names(source_cols)!='Other'],source_cols[names(source_cols)=='Other'])
  matrix_tracks<-list(data=data.frame(Source=source_value,stringsAsFactors=FALSE),
    cols=list(Source=source_cols),export_data=data.frame(Source=source_value,stringsAsFactors=FALSE))
  export(s,data.frame(matrix_order=seq_len(nrow(D)),isolate_id=rownames(D),
    run_accession=matrix_meta$run_accession,Source=source_value,check.names=FALSE),
    'Malawi_SNP_isolate_order_metadata')
  save_drawing(s,'Supplement_Malawi_pairwise_SNP_heatmap',
    snp_heatmap_draw(D,matrix_tracks,FALSE,max(D)),240,190,supp=TRUE)
  caption(s,'Figure_1A_Malawi_SNP_network',paste('Minimum-spanning network recalculated deterministically from the supplied pairwise SNP matrix.',
    if(s$cfg$collapse_zero)'Zero-distance groups were collapsed only after checking identical full distance profiles.' else 'Every isolate is a separate node, including zero-distance pairs.',
    'Node area reflects isolate count; sectors show the observed source fractions. Exact membership is exported for every node.',
    'Every edge is labelled numerically. The connected MST retains all nodes and differs in layout from the disconnected original illustration.',
    'Alternative equally minimal MSTs can exist; edge geometry and layout do not imply evolution or direction of transmission.'))
  note(s,'Figure 1A is a recomputed connected MST with the same input SNP distances, not a pixel replica of the original grouped/disconnected drawing. Review the new caption.')
}

module_map <- function(s) {
  need(c('sf','ggplot2'))
  file<-input_path(s,'map');stopif(!file.exists(file),paste('Missing',file,
    '- the map cannot be rebuilt from a manuscript raster image.'))
  layers<-sf::st_layers(file)$name;required<-c('boundary','waterways','sampling_sites','hospital')
  stopif(!all(required%in%layers),paste('GeoPackage layers required:',paste(required,collapse=', ')))
  objects<-setNames(lapply(required,function(k)sf::st_read(file,layer=k,quiet=TRUE)),required)
  for(k in required) {
    stopif(is.na(sf::st_crs(objects[[k]])),paste('Missing CRS:',k))
    stopif(any(!sf::st_is_valid(objects[[k]])),paste('Invalid geometries:',k))
    objects[[k]]<-sf::st_transform(objects[[k]],as.integer(s$cfg$map_epsg))
  }
  s$used<-unique(c(s$used,file));b<-objects$boundary;w<-objects$waterways;p<-objects$sampling_sites;h<-objects$hospital
  require_columns(p,c('site_id','category'),'Sampling sites');require_columns(w,'waterway_name','Waterways');require_columns(h,'label','Hospital')
  unique_ids(p$site_id,'Map site IDs')
  key<-identify(read_input(s,'site_members'))
  require_columns(key,'site_id','Site membership');stopif(any(!key$site_id%in%p$site_id),'An isolate maps to an unknown site')
  key<-attach_accessions(s,key);export(s,key,'Figure_1B_site_isolate_run_membership')
  bb<-sf::st_bbox(b);span<-as.numeric(bb['xmax']-bb['xmin'])
  scale_m<-pretty(c(0,span/4),n=3)[2];stopif(scale_m<=0,'Invalid map extent')
  x0<-as.numeric(bb['xmin'])+.06*span;y0<-as.numeric(bb['ymin'])+.05*as.numeric(bb['ymax']-bb['ymin'])
  p$category<-as_other(p$category);categories<-sort(unique(p$category));stopif(any(blank(categories)),'Map categories are missing')
  pal<-colours(categories);named<-w[!blank(w$waterway_name),,drop=FALSE]
  plot<-ggplot2::ggplot()+
    ggplot2::geom_sf(data=b,fill='#F3F6F2',colour='#B7C5BD',linewidth=.35)+
    ggplot2::geom_sf(data=w,colour='#81B0BF',linewidth=.35)+
    ggplot2::geom_sf(data=p,ggplot2::aes(colour=category),size=2.2)+
    ggplot2::geom_sf_text(data=p,ggplot2::aes(label=site_id),nudge_y=span*.009,size=2.1,colour=INK)+
    ggplot2::geom_sf(data=h,shape=8,size=3.0,colour='#AE5144')+
    ggplot2::geom_sf_text(data=h,ggplot2::aes(label=label),nudge_y=span*.024,size=2.5,colour=INK)+
    ggplot2::geom_sf_text(data=named,ggplot2::aes(label=waterway_name),size=2.0,colour='#517D8B',check_overlap=TRUE)+
    ggplot2::annotate('segment',x=x0,xend=x0+scale_m,y=y0,yend=y0,linewidth=.8)+
    ggplot2::annotate('text',x=x0+scale_m/2,y=y0+span*.016,label=paste(scale_m/1000,'km'),size=2.3)+
    ggplot2::scale_colour_manual(values=pal,name=NULL)+ggplot2::coord_sf(datum=NA)+
    ggplot2::theme_void()+pub_theme()+ggplot2::theme(axis.text=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),panel.grid=ggplot2::element_blank())
  save_gg(s,plot,'Figure_1B_Blantyre_map',180,165)
  caption(s,'Figure_1B_Blantyre_map','Map reconstructed from supplied vector GIS layers and explicit site-isolate membership. Site IDs, hospital labels and available waterway names are shown. Line widths do not encode an invented flow measurement. Exact source coordinates remain in the input GeoPackage; check disclosure permissions before public release.')
}

module_burden <- function(s) {
  need('ggplot2');b<-read_input(s,'burden')
  require_columns(b,c('isolate_id','group','n_amr_genes','snapshot_date'),'NCBI count comparison')
  unique_ids(b$isolate_id,'NCBI comparison IDs');b$n_amr_genes<-numeric_checked(b$n_amr_genes,'NCBI AMR counts')
  stopif(anyNA(b$n_amr_genes)||any(b$n_amr_genes<0|b$n_amr_genes!=round(b$n_amr_genes)),'NCBI gene counts must be known non-negative integers')
  stopif(any(blank(b$group))||length(unique(b$snapshot_date))!=1||any(blank(b$snapshot_date)),
    'Group and one documented snapshot_date are required')
  b$group<-as_other(b$group)
  sm<-do.call(rbind,lapply(split(b,b$group),function(x)data.frame(group=x$group[1],n=nrow(x),
    median=median(x$n_amr_genes),q25=unname(quantile(x$n_amr_genes,.25)),q75=unname(quantile(x$n_amr_genes,.75)),maximum=max(x$n_amr_genes))))
  order<-unique(b$group);b$group<-factor(b$group,levels=order);sm$group<-factor(sm$group,levels=order)

  p<-ggplot2::ggplot(b,ggplot2::aes(x=group,y=n_amr_genes,fill=group))+
    ggplot2::geom_violin(trim=TRUE,scale='width',alpha=.85,linewidth=.4,colour='white')+
    ggplot2::geom_linerange(data=sm,ggplot2::aes(x=group,y=median,ymin=q25,ymax=q75),inherit.aes=FALSE,linewidth=1.3,colour=INK)+
    ggplot2::geom_point(data=sm,ggplot2::aes(x=group,y=median),inherit.aes=FALSE,shape=21,size=2.5,fill='white',colour=INK)+
    ggplot2::geom_text(data=sm,ggplot2::aes(x=group,y=maximum,label=paste0('n = ',format(n,big.mark=','))),inherit.aes=FALSE,vjust=-1,size=2.8,colour=INK)+
    ggplot2::scale_fill_manual(values=colours(as.character(b$group)))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.02,.15)))+
    ggplot2::labs(x=NULL,y='Recorded AMR-gene count')+
    pub_theme()+ggplot2::theme(legend.position='none')
  save_gg(s,p,'Figure_3B_AMR_gene_counts',180,120)
  export(s,sm,'Figure_3B_summary');export(s,b,'Figure_3B_isolate_source_data')
  caption(s,'Figure_3B_AMR_gene_counts','Distribution of recorded AMR-gene counts in the supplied frozen NCBI comparison export. Denominators, medians and interquartile ranges are recomputed from that export, not copied from the manuscript or study-genome collection. The plot is descriptive and does not assign a clinical resistance category.')
}

module_plasmids <- function(s) {
  need(c('ggplot2','ggrepel'))
  p<-read_input(s,'plasmids');f<-read_input(s,'features');r<-read_input(s,'regions')
  require_columns(p,c('plasmid_id','isolate_id','length_bp','replicon','assembly_accession'),'Plasmid register')
  require_columns(f,c('plasmid_id','start','end','strand','gene','is_amr','region_id'),'Plasmid features')
  require_columns(r,c('plasmid_id','region_id','start','end'),'Resistance regions')
  unique_ids(p$plasmid_id,'Plasmid IDs');unique_ids(paste(r$plasmid_id,r$region_id,sep='\t'),'Plasmid-region IDs')
  p$length_bp<-numeric_checked(p$length_bp,'Plasmid length')
  stopif(anyNA(p$length_bp)||any(p$length_bp<=0|p$length_bp!=round(p$length_bp)),'Invalid plasmid lengths')
  stopif(any(blank(p$assembly_accession))||any(grepl('^[ESD]RR',p$assembly_accession)),
    'Supply closed-molecule sequence/assembly accessions, not read runs')
  for(k in c('start','end')) {f[[k]]<-numeric_checked(f[[k]],paste('feature',k));r[[k]]<-numeric_checked(r[[k]],paste('region',k))}
  check_coords<-function(x,what) {
    stopif(any(!x$plasmid_id%in%p$plasmid_id),paste('Unknown plasmid in',what))
    stopif(anyNA(x$start)||anyNA(x$end)||any(x$start<1|x$end<x$start|
      x$end>p$length_bp[match(x$plasmid_id,p$plasmid_id)]),paste('Invalid coordinates:',what))
  };check_coords(f,'features');check_coords(r,'regions')
  stopif(any(!f$strand%in%c('+','-')),'Gene strand must be + or -')
  f$is_amr<-numeric_checked(f$is_amr,'is_amr');stopif(anyNA(f$is_amr)||any(!f$is_amr%in%0:1),'is_amr must be explicit 0/1')
  p$label<-paste(p$isolate_id,p$plasmid_id,p$replicon,p$assembly_accession,sep=' | ')
  p$row<-rev(seq_len(nrow(p)));r$row<-p$row[match(r$plasmid_id,p$plasmid_id)]
  f$row<-p$row[match(f$plasmid_id,p$plasmid_id)];rc<-colours(r$region_id)
  top<-ggplot2::ggplot()+
    ggplot2::geom_segment(data=p,ggplot2::aes(x=0,xend=length_bp,y=row,yend=row),linewidth=2,colour='#D5DFE5')+
    ggplot2::geom_segment(data=r,ggplot2::aes(x=start,xend=end,y=row,yend=row,colour=region_id),linewidth=3.8)+
    ggplot2::scale_colour_manual(values=rc,name='Region')+
    ggplot2::scale_x_continuous(labels=function(z)z/1000,expand=ggplot2::expansion(mult=c(.01,.05)))+
    ggplot2::scale_y_continuous(breaks=p$row,labels=p$label)+
    ggplot2::labs(x='Position on the supplied closed plasmid (kb)',y=NULL)+
    pub_theme()+ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(size=7))
  if(has(s,'links')) {
    l<-read_input(s,'links');require_columns(l,c('query_plasmid','query_start','query_end','target_plasmid','target_start','target_end','percent_identity'),'Plasmid alignment blocks')
    for(k in c('query_start','query_end','target_start','target_end','percent_identity'))l[[k]]<-numeric_checked(l[[k]],k)
    stopif(anyNA(l)||any(l$percent_identity<0|l$percent_identity>100),'Invalid alignment block values')
    stopif(any(!l$query_plasmid%in%p$plasmid_id)|any(!l$target_plasmid%in%p$plasmid_id),'Unknown aligned plasmid')
    for(side in c('query','target')) {
      lengths<-p$length_bp[match(l[[paste0(side,'_plasmid')]],p$plasmid_id)]
      for(end in c('start','end')) {
        v<-l[[paste0(side,'_',end)]]
        stopif(any(v<1|v>lengths|v!=round(v)),'Alignment coordinates exceed a molecule or are not integer bp')
      }
    }
    poly<-do.call(rbind,lapply(seq_len(nrow(l)),function(i)data.frame(block=i,
      x=c(l$query_start[i],l$query_end[i],l$target_end[i],l$target_start[i]),
      y=c(rep(p$row[match(l$query_plasmid[i],p$plasmid_id)],2),rep(p$row[match(l$target_plasmid[i],p$plasmid_id)],2)))))

    top<-top+ggplot2::geom_polygon(data=poly,ggplot2::aes(x=x,y=y,group=block),fill='#7AA2B4',alpha=.13,colour=NA)
    export(s,l,'Figure_5_alignment_blocks')
  } else note(s,'No plasmid_alignments.csv: no sequence-homology ribbons are drawn. Backbone/detail plots do not reproduce an absent alignment.')
  save_gg(s,top,'Figure_5A_plasmid_backbones',240,max(110,35+nrow(p)*16))
  keys<-paste(r$plasmid_id,r$region_id,sep='\t');fk<-paste(f$plasmid_id,f$region_id,sep='\t')
  keep<-!blank(f$region_id);stopif(any(keep&!fk%in%keys),'A feature is assigned to an unknown resistance region')
  fi<-f[keep,,drop=FALSE];ri<-match(fk[keep],keys)
  stopif(any(fi$start<r$start[ri]|fi$end>r$end[ri]),'Feature falls outside its declared resistance region')
  stopif(!nrow(fi),'No resistance-region features to plot')
  fi$panel<-paste(fi$region_id,fi$plasmid_id,sep=' | ')
  fi$x0<-ifelse(fi$strand=='+',fi$start,fi$end);fi$x1<-ifelse(fi$strand=='+',fi$end,fi$start)
  fi$gene_key<-ifelse(blank(fi$gene),'Unnamed feature',fi$gene)
  pal<-colours(fi$gene_key)
  detail<-ggplot2::ggplot(fi,ggplot2::aes(y=0))+
    ggplot2::geom_segment(ggplot2::aes(x=x0,xend=x1,yend=0,colour=gene_key),linewidth=1.4,
      arrow=grid::arrow(type='closed',length=grid::unit(1.5,'mm')))+
    ggrepel::geom_text_repel(ggplot2::aes(x=(start+end)/2,label=gene_key),direction='y',nudge_y=.45,
      min.segment.length=0,max.overlaps=Inf,max.time=Inf,max.iter=20000,size=2.4,seed=s$cfg$seed)+
    ggplot2::facet_wrap(~panel,scales='free_x',ncol=2)+ggplot2::scale_colour_manual(values=pal)+
    ggplot2::scale_x_continuous(labels=function(z)z/1000)+ggplot2::coord_cartesian(ylim=c(-.2,1.6),clip='off')+
    ggplot2::labs(x='Position on the supplied plasmid (kb)',y=NULL)+
    pub_theme()+ggplot2::theme(legend.position='none',axis.text.y=ggplot2::element_blank(),panel.grid=ggplot2::element_blank())
  save_gg(s,detail,'Figure_5B_resistance_region_details',240,max(140,ceiling(length(unique(fi$panel))/2)*52))
  export(s,p,'Figure_5_plasmid_register');export(s,f,'Figure_5_feature_coordinates');export(s,r,'Figure_5_region_coordinates')
  caption(s,'Figure_5',paste('Resolved plasmid backbones with isolate/molecule accession labels, followed by enlarged annotated resistance regions.',
    'Coordinates, orientations, region assignments and lengths come only from supplied tables. Identical gene names use the same colour.',
    'Homology polygons are included only when an explicit alignment-block table is supplied.',
    'These structural comparisons do not establish a particular cointegration or transfer event.'))
}

self_tests <- function() {
  stopifnot(identical(clean_id(c('ABC_assembly_contigs.fasta','CIV13XE8_S20_L001','44341_3#96')),
    c('ABC','CIV13XE8','44341_3_96')))
  d<-matrix(c(0,2,2,0),2,dimnames=list(c('a','b'),c('a','b')))
  stopifnot(identical(validate_distance(d[,2:1,drop=FALSE]),d))
  stopifnot(inherits(try(validate_distance(matrix(c(0,2,3,0),2,dimnames=dimnames(d))),silent=TRUE),'try-error'))
  stopifnot(inherits(try(unique_ids(c('a','a'),'test'),silent=TRUE),'try-error'))
  z<-matrix(c(0,0,2,0,0,2,2,2,0),3,dimnames=list(c('a','b','c'),c('a','b','c')))
  g<-collapse_zero(z);stopifnot(nrow(g$D)==2,nrow(g$membership)==3,nrow(kruskal(g$D))==1)
  stopifnot(is.na(numeric_checked(NA_character_,'test')))
  stopifnot(identical(as_other(c(NA,'','Unknown','Malawi')),c('Other','Other','Other','Malawi')))
  stopifnot(identical(canonical_runs('SRR2;SRR1;SRR1'),'SRR1;SRR2'))
  message('PASS: ID, missing-value, SNP-matrix, accession and zero-distance grouping self-tests.')
  invisible(TRUE)
}
finalise <- function(s) {
  utils::write.csv(s$status,file.path(s$out,'audit','build_status.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(s$out,'audit','sessionInfo.txt'))
  actual<-s$used[file.exists(s$used)]
  if(length(actual)) {
    hash<-data.frame(path=actual,md5=unname(tools::md5sum(actual)))
    if(requireNamespace('digest',quietly=TRUE))hash$sha256<-vapply(actual,function(x)digest::digest(file=x,algo='sha256'),character(1))
    utils::write.csv(hash,file.path(s$out,'audit','input_checksums.csv'),row.names=FALSE)
  }
  products<-list.files(s$out,recursive=TRUE,full.names=TRUE)
  if(length(products))utils::write.csv(data.frame(file=substring(products,nchar(s$out)+2),bytes=file.info(products)$size,
    md5=unname(tools::md5sum(products))),file.path(s$out,'audit','output_checksums.csv'),row.names=FALSE)
  blocked<-any(s$status$state=='BLOCKED')
  label<-if(blocked)'PARTIAL BUILD - inspect build_status.csv' else 'SELECTED PANELS BUILT - scientific review still required'
  writeLines(c(label,'A successful render is not laboratory, biological, public-access or manuscript-validation approval.'),
    file.path(s$out,'BUILD_STATUS.txt'))
  message('\n',label,'\nOutputs: ',s$out)
  invisible(s)
}
