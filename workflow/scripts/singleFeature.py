def single_feature(seqid: str, end: int, start: int = 1, source: str ="user") -> list:
    '''
    Generate a single transcript/exon gene feature with predefined attributes.
    '''

    import gffutils
    gene = gffutils.Feature(seqid=seqid,
                            source=source,
                            featuretype='gene',
                            start=start,
                            end=end,
                            strand='+',
                            id=f"{seqid}_gene",
                            attributes={'gene_id': [f"{seqid}_gene"], 
                                        'gene_type': ["synthetic_construct"],
                                        'gene_name': [f"{seqid}_gene"],
                                        'level': ["2"],
                                        'tag': ["synthetic"]})
    transcript = gffutils.Feature(seqid=seqid,
                            source=source,
                            featuretype='transcript',
                            start=start,
                            end=end,
                            strand='+',
                            id=f"{seqid}_transcript",
                            attributes={'gene_id': [f"{seqid}_gene"],
                                        'transcript_id': [f"{seqid}_transcript"],
                                        'gene_type': ["synthetic_construct"],
                                        'gene_name': [f"{seqid}_gene"],
                                        'transcript_type': ["synthetic"],
                                        'transcript_name': [f"{seqid}_transcript"],
                                        'level': ["2"],
                                        'tag': ["synthetic"],
                                        'Parent': [f"{seqid}_gene"]})
    exon = gffutils.Feature(seqid=seqid,
                            source=source,
                            featuretype='exon',
                            start=start,
                            end=end,
                            strand='+',
                            id=f"{seqid}_exon",
                            attributes={'gene_id': [f"{seqid}_gene"], 
                                        'transcript_id': [f"{seqid}_transcript"],
                                        'gene_type': ["synthetic_construct"],
                                        'gene_name': [f"{seqid}_gene"],
                                        'transcript_type': ["synthetic"],
                                        'transcript_name': [f"{seqid}_transcript"],
                                        'exon_number': ["1"],
                                        'exon_id': [f"{seqid}_exon"],
                                        'level': ["2"],
                                        'tag': ["synthetic"],
                                        'Parent': [f"{seqid}_transcript"]})
    return [gene, transcript, exon]

def format_attributes(attrs):
    formatted = []
    for key, values in attrs.items():
        for value in values:
            formatted.append(f'{key} {value}') if key in ["level", "exon_number"] else formatted.append(f'{key} "{value}"')
    return '; '.join(formatted) + ';'