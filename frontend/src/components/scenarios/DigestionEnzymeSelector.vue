<template lang="pug">

div
  h1 Digestion Selector
  img.icon.center-block(slot='title-img', :src='infos.icon' title='NOT a proto-nazi symbol !')
  p.center.
    Find the best enzymes to digest your constructs, for verification or
    identification purposes.

  p.
    If no single digestion works for all constructs,
    2 or 3 digestions that collectively cover all constructs will be suggested.

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
    sequencesuploader(v-model='form.files')

    h4.formlabel Ladder
    el-select(v-model='form.ladder', placeholder='Select')
      el-option(v-for='item in ladder_options', :label='item.label',
                :value='item.value', :key='item.value')

    p.inline Bands precision (%)
      el-input-number.inline(v-model="form.bandsPrecision", size="small",
                             :min=1, :max=100)
      helper(help='Precision in (%) of the total migration distance span. Determines how far away bands should be to get distinguished.')

    h4.formlabel Possible enzymes

    p.loadData(@click='form.possibleEnzymes = enzymesPreselection') &#9656; Click me for ~60 preset enzymes.
    el-input(type='textarea', :rows='6', placeholder='Example: EcoRI, XbaI, XhoI, ...',
             v-model="form.possibleEnzymes")

    p.inline Max. enzymes in one digestion
      el-input-number.inline(v-model="form.maxEnzymes", size="small", :min='1', :max='3')

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

</template>

<script>
import learnmore from '../../components/widgets/LearnMore'
import sequencesuploader from '../../components/widgets/SequencesUploader'
import digestionset from '../../components/widgets/DigestionSelectorSet'

var infos = {
  title: 'Select Digestions',
  navbarTitle: 'Select digestions',
  path: 'digestion-selector',
  description: '',
  backendUrl: 'start/select_digestions',
  icon: require('assets/images/select_digestion.svg')
}

export default {
  data: function () {
    return {
      form: {
        ladder: '100_to_4k',
        possibleEnzymes: 'EcoRI, XbaI, XhoI',
        maxEnzymes: 1,
        bandsPrecision: 5,
        make_report: false,
        goal: 'ideal',
        bandsRange: [3, 5],
        files: []
      },
      infos: infos,
      ladder_options: [
        {
          label: 'Ladder 100 bp - 4000 bp',
          value: '100_to_4k'
        },
        {
          label: 'Ladder 35 bp - 5 kbp (AATI)',
          value: '35_to_5k'
        },
        {
          label: 'Ladder 75 bp - 15 kbp (AATI)',
          value: '75_to_15k'
        }
      ],
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
        'AatII, Acc65I, AccI, AflII, AluI, ApaI, ApaLI, AseI, AvaI, AvaII, ' +
        'BamHI, BanI, BanII, BclI, BglII, BsmI, Bsp1286I, BssHII, ClaI, DdeI, ' +
        'DpnI, DraI, EcoRI, EcoRV, FokI, HaeII, HaeIII, HhaI, HincII, ' +
        'HindIII, HinfI, HpaI, HpaII, KpnI, MboI, MluI, MnlI, MspA1I, MspI, ' +
        'NaeI, NarI, NciI, NcoI, NdeI, NheI, NruI, NsiI, PspOMI, PstI, PvuI, ' +
        'PvuII, RsaI, SacI, SacII, SalI, Sau96I, SmaI, SnaBI, SpeI, SphI, ' +
        'SspI, StuI, StyI, XbaI, XhoI, XmaI'
      )
    }
  },
  components: {
    sequencesuploader,
    learnmore,
    digestionset
  },
  infos: infos,
  methods: {
    handleSuccess: function (evt) {
      console.log(evt)
    },
    validateForm: function () {
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
