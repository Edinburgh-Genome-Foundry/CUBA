<template lang="pug">
.page
  h1 {{infos.title}}
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Get enzyme suggestions for your restriction digests:",
            :tweetUrl="'https://cuba.genomefoundry.org/' + infos.path")
  img.icon.center-block(slot='title-img', :src='infos.icon')
  p.scenario-description.
    Find optimal Sanger sequencing plans for your constructs, maximizing primer
    reuse and minimizing the total numbers of reads and of new primers to buy.

  .form
    h4.formlabel Validation type
    el-select(v-model='form.goal', placeholder='Select a goal')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')

    h4.formlabel Strands covered
    el-select(v-model='form.strand', placeholder='Select a strand')
      el-option(v-for='item in strand_options', :label='item.label',
                :value='item.value', :key='item.value')

    h4.formlabel Constructs Sequences
      helper.
        The sequences should be Genbanks with features annotated "cover" for the
        zones to be covered (it can be one big zone covering the full sequence)
        or "no_primer" to forbid primers in certain zones.
    collapsible(title='Examples')
      file-example(filename='20_random_constructs.zip',
                   @input='function (e) {form.constructs.push(e)}',
                   fileHref='/static/file_examples/select_primers/20_random_constructs.zip',
                   imgSrc='/static/file_examples/select_primers/20_random_constructs.png')
        p.
          Collection of constructs assembled using random parts from the EMMA standard.
          Unzip the file and drag the genbank files into the file upload area.

    files-uploader(v-model='form.constructs', tip='Fasta or Genbank files')
    el-checkbox(v-model='form.circularSequences') Sequences are circular

    h4.formlabel Available primers
    collapsible(title='Examples')
      file-example(filename='available_primers.fa',
                   @input='function (e) {form.availablePrimers.push(e)}',
                   fileHref='/static/file_examples/select_primers/available_primers.fa',
                   imgSrc='/static/file_examples/select_primers/available_primers.png')
        p.
          Collection of constructs assembled using random parts from the EMMA standard.
          Unzip the file and drag the genbank files into the file upload area.

    files-uploader(v-model='form.availablePrimers', help='Fasta or Genbank files')

    h4.formlabel Primers properties

    p.inline Ideal read range: from <b>{{ form.readRange[0] }}bp</b> after primer annealing to <b>{{ form.readRange[1] }}bp</b>
      el-slider(range v-if="" v-model='form.readRange',
                :min='0', :max='2000', :step='50')
    p.inline Primer melting temperatures from <b>{{ form.tmRange[0] }}C</b> to <b>{{ form.tmRange[1] }}C</b>
      el-slider(range v-if="" v-model='form.tmRange',
                :min='20', :max='100', :step='1')
    p.inline Digits in new primers names:
      el-input-number.inline(v-model="form.newPrimerDigits", size="small",
                             :min='1', :max='6')

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Select primers',
                    v-model='queryStatus')
    progress-bars(:bars='queryStatus.polling.data.bars', :order="['record', 'primer']"
                  v-if='queryStatus.polling.inProgress && queryStatus.polling.data')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")

    .results(v-if='!queryStatus.polling.inProgress')
      p.results-summary(v-if='queryStatus.result.summary',
                        v-html="queryStatus.result.summary")
      download-button(v-if='queryStatus.result.zip_file',
                      :filedata='queryStatus.result.zip_file')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>

var infos = {
  title: 'Select Primers',
  navbarTitle: 'Select primers',
  path: 'select_primers',
  description: '',
  backendUrl: 'start/select_primers',
  icon: require('../../assets/images/select_primers.svg'),
  poweredby: ['primavera']
}

export default {
  data () {
    return {
      form: {
        goal: 'sanger_sequencing',
        circularSequences: true,
        constructs: [],
        readRange: [150, 800],
        tmRange: [55, 70],
        strand: 'both',
        newPrimerDigits: 4,
        availablePrimers: []
      },
      infos: infos,
      goal_options: [
        {
          label: 'Sanger Sequencing',
          value: 'sanger_sequencing'
        },
        {
          label: 'Validation PCR (NOT IMPLEMENTED YET)',
          value: 'validation_pcr'
        }
      ],
      strand_options: [
        {
          label: 'Both',
          value: 'both'
        },
        {
          label: '3\'-5\'',
          value: '3-5'
        },
        {
          label: '5\'-3\'',
          value: '5-3'
        },
        {
          label: 'Any',
          value: 'any'
        }
      ],
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      }
    }
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (this.form.constructs.length < 1) {
        errors.push('Provide at least one construct !')
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

.el-input-number.inline {
  margin-bottom: -9px;
  margin-left: 10px;
  width: 130px;
}

p.loadData {
  font-size: 0.8em;
  margin-bottom: 0;
  cursor: pointer;
}

.bands-range {
  .el-slider {
    display: inline-block;
    width: 150px;
    margin-bottom: -12px;
    margin-left: 10px;
  }
}
</style>
