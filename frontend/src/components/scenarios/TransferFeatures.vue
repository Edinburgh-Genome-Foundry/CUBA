<template lang="pug">

.page
  h1  {{ infos.title }}
  img.icon.center-block(slot='title-img', :src='infos.icon')
  
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Easily transfer features between Genbanks:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path"
            style='margin-top: 2em;')
  .scenario-description
    p.
      Submit several genbank records, get new versions where the features
      of each record are copied to every records on locations with identical 
      subsequence.
    p.
      This is useful to transfer features from a fully-annotated record to an
      un-annotated record
      
    p.
      For instance, many DNA manufacturers provide unannotated
      plasmid records when you order genetic parts. This app enables to
      manufacturer records using the annotations from your original part records.

  

  

  .form

    h4.formlabel Genbank records
    collapsible(title='Examples')
      file-example(filename='Part and un-annotated plasmid',
                   fileHref='/static/file_examples/transfer_features/part_and_unannotated_plasmid.zip',
                   @input='function (e) { form.files.push(e) }'
                   imgSrc='/static/file_examples/generic_logos/sequences_records.png')
        p
          :markdown-it
            This ZIP file contains two genbank records:

            - A record of an annotated part, containing an expression module.
            - A record of a plasmid which contains the part but the part was not properly annotated

            The app will automatically detect where the part is located in
            the plasmid and automatically copy the features from the part record to the
            plasmid record.

    files-uploader(v-model='form.files', :multiple='true',
                      text="Drop multiple Genbank/Fasta (or click to select)")
    p.inline
      el-checkbox(v-model='form.use_filenames_as_ids') Use file names as record IDs
    p.inline Minimal homology block size (nucleotides):
      el-input-number.inline(v-model="form.min_block_size", size="small", :min=10, :max=5000)

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Transfer features',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress && queryStatus.result')
      download-button(v-if='queryStatus.result.zip_file',
                      :filedata='queryStatus.result.zip_file',
                      text='Download Sequences')
      ul
        li(v-for="[name, n] in queryStatus.result.added_features", :key='name').
          #[b {{n}}] features added to record #[b {{name}}]

  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'

var infos = {
  title: 'Transfer Genbank Features',
  navbarTitle: 'Transfer Genbank Features',
  path: 'transfer-features',
  description: '',
  backendUrl: 'start/transfer_features',
  icon: require('../../assets/images/transfer_features.svg'),
  poweredby: ['geneblocks']
}

export default {
  data () {
    return {
      form: {
        files: [],
        min_block_size: 50,
        use_filenames_as_ids: true
      },
      infos: infos,
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  components: {
    learnmore
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (this.form.files.length === 0) {
        errors.push('Provide at least two sequences.')
      }
      return errors
    }
  }
}
</script>

<style lang='scss' scoped>

h4.formlabel {
  text-align: center;
  text-transform: uppercase;
  margin-top: 40px
}

.form {
  margin: 50px auto;
  max-width: 500px;
}

.title-img {
  height:80px;
  margin-top: -20px;
  margin-bottom: 20px;

}

.el-checkbox {
  font-weight: normal;
}

.el-select {
  width: 100%
}

.report-radio {
  text-align: center;
  margin: 0 auto;
  .el-radio {
    margin-bottom: 20px;
  }
}

.files-number p {
  margin-top: -10px;
  font-size: 0.8em;
}

.figure-preview {
  margin-bottom: 5em;
  img {
    width:100%;
  }
}

.scenario-description {
  max-width: 800px;
}
</style>
