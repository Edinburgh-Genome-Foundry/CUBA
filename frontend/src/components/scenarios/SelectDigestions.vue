<template lang="pug">
.page
  h1 {{infos.title}}
  img.icon.center-block(slot='title-img', :src='infos.icon' title='NOT a proto-nazi symbol !')
  p.center.
    Find the best enzymes to digest your constructs, for verification or
    identification purposes.
  p.
    If no single digestion works for all constructs,
    2 or 3 digestions that collectively cover all constructs will be suggested.
  web-links(:emailSubject="'[CUBA] Feedback on web app: ' + infos.title",
            tweetMessage="Get enzyme suggestions for your restriction digests:",
            :tweetUrl="'http://cuba.genomefoundry.org/' + infos.path")



  .form
    h4.formlabel The digestions should produce
    el-select(v-model='form.goal', placeholder='Select')
      el-option(v-for='item in goal_options', :label='item.label',
                :value='item.value', :key='item.value')

    .bands-range(v-if="form.goal === 'ideal'")
      p.inline Bands in <i>sensitivity</i> zone: {{ form.bandsRange[0] }} to {{ form.bandsRange[1] }}
        el-slider(range show-stops v-if="" v-model='form.bandsRange',
                       :max=10, :min=1)
        helper(help='Bands in the high-band-size-sensitivity zone, between 10% and 90% of the maximal migration distance.')

    h4.formlabel Constructs Sequences


    files-uploader(v-model='form.files', help='Fasta or Genbank files')
    el-checkbox(v-model='form.circularSequences') Sequences are circular

    h4.formlabel Ladder
    ladderselector(v-model='form.ladder')

    h4.formlabel Possible enzymes

    p.loadData(@click='form.possibleEnzymes = enzymesPreselection') &#9656; Click me for ~60 preset enzymes.
    el-input(type='textarea', :rows='6', placeholder='Example: EcoRI, XbaI, XhoI, ...',
             v-model="form.possibleEnzymes")

    p.inline Max. enzymes in one digestion
      el-input-number.inline(v-model="form.maxEnzymes", size="small", :min='1', :max='3')
    p.inline Number of digestions: &nbsp; &nbsp;
      el-radio(v-model='form.maxDigestions', label='1') 1
      el-radio(v-model='form.maxDigestions', label='2') 2

    el-checkbox(v-model='form.showBandsSizes') Show bands sizes in plot.

    backend-querier(:form='form', :backendUrl='infos.backendUrl',
                    :validateForm='validateForm', submitButtonText='Select digestions',
                    v-model='queryStatus')
    el-alert(v-if='queryStatus.requestError', :title="queryStatus.requestError",
       type="error", :closable="false")
    .results(v-if='!queryStatus.polling.inProgress && queryStatus.result.digestions')
      .results-summary(v-if='queryStatus.result.digestions')
        p <b>Results</b>
        p(v-if='queryStatus.result.digestions.length').
          The solution found involves {{ queryStatus.result.digestions.length}} digestions
        p(v-else).
          No solution found :'(
      .center
        img(v-if='queryStatus.result.figure_data', :src='queryStatus.result.figure_data')
  powered-by(:softwareNames='infos.poweredby')
</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'
import ladderselector from '../../components/widgets/LadderSelector'

var infos = {
  title: 'Select Digestions',
  navbarTitle: 'Select digestions',
  path: 'select_digestions',
  description: '',
  backendUrl: 'start/select_digestions',
  icon: require('assets/images/select_digestion.svg'),
  poweredby: ['bandwitch', 'bandwagon']
}

export default {
  data () {
    return {
      form: {
        ladder: '100_to_4k',
        possibleEnzymes: 'EcoRI, XbaI, XhoI',
        maxEnzymes: 1,
        maxDigestions: '1',
        make_report: false,
        goal: 'ideal',
        bandsRange: [3, 5],
        circularSequences: true,
        showBandsSizes: false,
        files: []
      },
      infos: infos,
      goal_options: [
        {
          label: 'Good patterns for all constructs',
          value: 'ideal'
        },
        {
          label: 'Different patterns for all constructs, for identification',
          value: 'separating'
        }
      ],
      queryStatus: {
        polling: {},
        result: {},
        requestError: ''
      },
      enzymesPreselection: (
        'AccI, AclI, AflII, AflIII, AgeI, ApaLI, AvaI, BlnI, BmtI, BsmI, ' +
        'BssHII, DdeI, DraI, Eco47III, EcoRI, EcoRV, HindII, HindIII, HinfI, ' +
        'HpaI, MluI, MspA1I, MunI, NaeI, NciI, NcoI, NdeI, NsiI, PstI, PvuII, ' +
        'SacI, SacII, SalI, ScaI, SfaNI, SnaBI, SpeI, SphI, SspI, StyI, VspI, ' +
        'XhoI, XmaI, ZraI'
      )
    }
  },
  components: {
    sequencesuploader,
    learnmore,
    digestionset,
    ladderselector
  },
  infos: infos,
  methods: {
    handleSuccess (evt) {
      console.log(evt)
    },
    validateForm () {
      var errors = []
      if (this.form.possibleEnzymes.length < 2) {
        errors.push('Provide at least two different restriction enzymes.')
      }
      if (this.form.files.possibleEnzymes === 0) {
        errors.push('Provide at least one construct file')
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
