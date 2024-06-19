
     ##new.concats.dt = gr2dt(new.concats)
library(devtools)
library(gUtils)
library(data.table)
devtools::load_all('~/git/chromunity')

args = commandArgs(trailingOnly = TRUE)

concatemers = readRDS(args[1])

concatemers$cid = concatemers$read_idx

if (class(concatemers)[[1]] == 'GRanges') {
    new.concats.dt = gr2dt(concatemers)
    concatemers.gr = concatemers
} else {
    colnames(concatemers)[colnames(concatemers) == 'chrom'] = 'seqnames'
    concatemers.gr = dt2gr(concatemers)
    new.concats.dt = concatemers
}

new.concats.dt[, bincount := .N, by='read_idx']

unique.concats = unique(new.concats.dt, by='read_idx')

    ##unique.concats = unique.rerun
    
    #number of chromosomes a concatemer hits    
new.concats.dt[, num.chr := length(unique(seqnames)), by='read_idx']

##number of times a concatemer hits the given chromosome
new.concats.dt[, chr.count := .N, by=c('read_idx','seqnames')]

##unique.chr.concats
unique.chr.concats = unique(new.concats.dt, by=c('read_idx','seqnames'))
unique.chr.concats[, intra.chr.contacts := choose(chr.count, 2)]

##unique.concats[, contact.count := choose(bincount, 2), by='read_name']
unique.chr.concats[, contacts.count := choose(bincount, 2), by='read_idx']
unique.concats = unique(unique.chr.concats, by=c('read_idx'))

intra.vs.inter = unique.chr.concats[, .(intra.contacts = sum(intra.chr.contacts), contacts.count), by='read_idx']
intra.vs.inter = unique(intra.vs.inter, by='read_idx')

total.inter.contacts = sum(intra.vs.inter$contacts.count) - sum(intra.vs.inter$intra.contacts)
intra.vs.inter = unique(intra.vs.inter, by='read_idx')

percent.cis = sum(intra.vs.inter$intra.contacts) / sum(intra.vs.inter$contacts.count)
##pairs.new.row$percent.cis = percent.cis

##contact count
contact.count = sum(unique.concats$contacts.count)
##pairs.new.row$contact.count = contact.count

##contacts per gb unique.c
##contacts.per.gb = contact.count / ((total.gb) / (10^9))
##pairs.new.row$contacts.per.gb = contacts.per.gb

##contacts.per.gb = (pairs.table[1]$contact.count) / (pairs.table[1]$total_bases / (10^9))
##contacts.per.gb
##pairs.table[1]$contacts.per.gb = contacts.per.gb

###contact order distribution
concat.analysis = data.table()
concat.analysis$lessone = sum(unique.concats$bincount == 1)

total.reads = dim(unique.concats)[[1]]

concat.analysis$multiway = sum(unique.concats$bincount > 1)
concat.analysis$two = sum(unique.concats$bincount == 2) / total.reads
concat.analysis$three = sum(unique.concats$bincount == 3)  / total.reads
concat.analysis$four = sum(unique.concats$bincount == 4) / total.reads
concat.analysis$fivetosix = sum((unique.concats$bincount == 5) | (unique.concats$bincount == 6)) / total.reads
concat.analysis$seventoeleven = sum((unique.concats$bincount >= 7) & (unique.concats$bincount <= 11)) / total.reads
concat.analysis$twelvetotwentyone = sum((unique.concats$bincount >= 12) & (unique.concats$bincount <= 21)) / total.reads
concat.analysis$twentytwotofortynine = sum((unique.concats$bincount >= 22) & (unique.concats$bincount <= 49)) / total.reads
concat.analysis$greaterthan50 = sum(unique.concats$bincount >= 50) / total.reads

concat.analysis$contact.count = contact.count
concat.analysis$percent.cis = percent.cis

###reporting binned concatemer contacts
###to disentangle dupage

print(concatemers.gr)
bins = gr.tile(hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38::Hsapiens"), width=10000)
concatemers.gr$binid = gr.match(concatemers.gr, bins, max.slice = 1e6, mc.cores = 5, verbose = FALSE)

concatemers.gr = concatemers.gr %Q% (!is.na(binid))

                               
concatemers.gr$cid = concatemers.gr$read_idx
concats.dt = as.data.table(concatemers.gr)[, `:=`(count, .N), by = cid]
concats.dt[, cidi := as.integer(cid)]
##dedupe bin hits

concats.dt = unique(concats.dt, by=c('cidi','binid'))
concats.dt[, bincount := .N, by='read_idx']
concats.dt[, contacts.count := choose(bincount, 2), by='read_idx']
deduped.contact.count = sum(unique(concats.dt, by='cidi')$contacts.count)
concat.analysis$deduped.contact.count = deduped.contact.count

saveRDS(concat.analysis, 'concatemer_analysis.rds')
