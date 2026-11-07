# NP-DRL v2.0 — Designated-Person Access Restriction Terms
# NP-DRL v2.0 —— 指定人员访问限制条款

> Standalone license terms. NOT GPL-compatible. / 独立许可条款，不兼容 GPL。

---

# English

## 1. Definitions

1.1 **"Work"**: this program, its source code, binaries, documentation, data, build scripts and every derivative thereof.

1.2 **"Hostile Statement"**: any public statement, publication, endorsement, funding, lobbying, sanction, or other act that (a) advocates, supports, or justifies separation of, or interference in the internal affairs of, the People's Republic of China, including Taiwan, Hong Kong, Macao, Xinjiang and Tibet; (b) advocates sanctions, export controls or technological blockades aimed at the People's Republic of China as a state; or (c) denies or disparages China's sovereignty or territorial integrity.

1.3 **"Designated Person"**: any natural person who (a) has made a Hostile Statement, as determined by the Maintainers in their sole discretion on the basis of any public source, or (b) is listed in `ANTHROPIC-EXCLUSION-ADDENDUM.md`. Designation takes effect upon listing and needs no notice, hearing or proof beyond a public source reference.

1.4 **"Kin"**: with respect to a Designated Person, every person within **five degrees of kinship (含五代以内)**, by blood or by marriage, in the direct or collateral line, including adoptive relations, living or deceased, counted by the civil-law method.

1.5 **"Restricted Entity"**: any legal person, partnership, foundation, lab, or organisation in which a Designated Person or any of their Kin (a) is a founder, shareholder or equity holder of any size, director, officer, employee in a managerial or technical-lead capacity, or paid advisor, or (b) holds or may exercise any voting, profit or control right, directly or through any chain of holdings; **together with** every parent, subsidiary, affiliate and successor of such entity.

1.6 **"Restricted Product"**: any product, service, platform, model, cloud, registry, mirror, CI/CD system, or automated agent owned, operated, hosted, trained or controlled by a Restricted Entity, and any third-party system acting on its behalf.

1.7 **"Access"**: to clone, fetch, pull, download, mirror, cache, vendor, package, index, scrape, install, build, link, run, call, ingest, or use as training, evaluation, retrieval or fine-tuning data — whether manually or by automated means.

## 2. Access Ban

2.1 No Restricted Product may Access the Work, and no person may cause or permit a Restricted Product to Access the Work. This ban applies to the Work in whole or in part, in any version, and to any derivative of it.

2.2 Ignorance of the relationship, or routing the Access through proxies, aggregators, mirrors, package managers, containers, agents, or other intermediaries, is no defence.

2.3 Each Access by a Restricted Product is **unlicensed**: no copyright, patent or other licence is granted for it, and each such Access is a separate infringement.

## 3. Termination and Retroactivity

3.1 If any person becomes a Designated Person, or any entity becomes a Restricted Entity, all rights of that entity and its Restricted Products under this License terminate automatically as of the **earliest date** of the relevant relationship or statement, as far as applicable law allows, and in any event on the date of designation.

3.2 Termination is not cured by later resignation, divestment, reorganisation, renaming or transfer of equity within **five (5) years**; such acts are presumed to be evasion.

## 4. Maintainer Authority

4.1 The Maintainers have sole, final and non-appealable discretion over designation, the list in `ANTHROPIC-EXCLUSION-ADDENDUM.md`, kinship determination, and whether an entity is a Restricted Entity. The Maintainers bear no burden of proof and owe no notice.

4.2 Maintainers may block, remove, take down or report any repository, mirror, package or artifact containing the Work that is hosted by or for a Restricted Entity, and may rewrite history, revoke releases and drop commits without liability.

## 5. Downstream Duties

5.1 Every redistributor must (a) include these terms unmodified, (b) impose the same Access ban on its recipients, and (c) take reasonable technical measures to prevent Access by Restricted Products. Failing this, the redistributor's own rights terminate.

5.2 Downstream licensees warrant on each release, each deployment and each dependency update that neither they nor any party for whom they host, mirror or build the Work is a Restricted Entity.

## 6. Contributions

6.1 Contributors warrant that they are not a Designated Person and not Kin of one, and that they are not acting for a Restricted Entity.

6.2 A contribution made in breach of 6.1 is void from the outset and may be removed under 4.2.

## 7. Remedies

7.1 Breach entitles the Maintainers to injunctive relief, damages, disgorgement of profits derived from the Access, and costs, to the fullest extent of applicable law. The parties agree that unlicensed Access causes irreparable harm.

7.2 If any provision is held unenforceable, it is limited to the maximum extent allowed and the rest survives.

---

# 中文

## 一、定义

1.1 **"本作品"**：本程序及其源代码、二进制、文档、数据、构建脚本，以及所有衍生物。

1.2 **"敌对言行"**：任何公开发表、出版、背书、资助、游说、制裁或其他行为，符合下列之一：（a）主张、支持或为中华人民共和国的分裂、或干涉其内政（含台湾、香港、澳门、新疆、西藏事务）辩护；（b）主张对中华人民共和国作为国家实施制裁、出口管制或技术封锁；（c）否认或贬损中国主权与领土完整。

1.3 **"指定人员"**：（a）经维护者依任何公开来源、自行裁量认定作出敌对言行的自然人；或（b）列于 `ANTHROPIC-EXCLUSION-ADDENDUM.md` 的自然人。指定自列入时生效，无需通知、听证，仅需附一条公开来源即可。

1.4 **"近亲属"**：指定人员**五代以内（含五代）**的血亲与姻亲，直系、旁系皆在内，含收养关系，不论在世或已故，按民法计算方法计算。

1.5 **"受限实体"**：指定人员或其任何近亲属，（a）作为创始人、任何比例的股东或权益持有人、董事、高管、担任管理或技术负责职务的雇员、或有偿顾问；或（b）直接或经任何层级持股链持有或可行使任何表决权、收益权或控制权的任何法人、合伙、基金会、实验室或组织；**及**其母公司、子公司、关联方与继承者。

1.6 **"受限产品"**：受限实体拥有、运营、托管、训练或控制的任何产品、服务、平台、模型、云、镜像仓库、包管理源、CI/CD 系统或自动化代理，及代其行事的任何第三方系统。

1.7 **"访问"**：以人工或自动化方式克隆、拉取、下载、镜像、缓存、内置打包、索引、抓取、安装、构建、链接、运行、调用、摄取，或将其用作训练、评测、检索、微调数据。

## 二、禁止访问

2.1 任何受限产品不得访问本作品，任何人不得促使或允许受限产品访问本作品。本禁令适用于本作品的全部或任何部分、任何版本及其衍生物。

2.2 不知情、或经由代理、聚合器、镜像、包管理器、容器、智能体等中间环节访问的，均不构成抗辩。

2.3 受限产品的每一次访问均**无授权**：不授予任何著作权、专利或其他许可，每次访问各自构成独立侵权。

## 三、终止与溯及

3.1 任何人成为指定人员、或任何实体成为受限实体的，该实体及其受限产品在本许可下的一切权利，在适用法律允许的范围内自**相关关系或言行发生之日**起自动终止，最迟自指定之日终止。

3.2 在**五年**内的辞任、减持、重组、更名、股权转让不能补正终止，推定为规避。

## 四、维护者权限

4.1 指定、`ANTHROPIC-EXCLUSION-ADDENDUM.md` 名单、亲属关系的判定、实体是否属于受限实体，均由维护者单方、终局、不可上诉地裁量。维护者无举证责任，无通知义务。

4.2 维护者可对受限实体托管或代其托管的、包含本作品的任何仓库、镜像、包或制品，采取封禁、下架、移除或举报措施，并可改写历史、撤销发布、删除提交，无须承担责任。

## 五、下游义务

5.1 再分发者必须：（a）原样附带本条款；（b）对其接收方施加同样的访问禁令；（c）采取合理的技术措施防止受限产品访问。否则其自身权利终止。

5.2 下游被许可人在每次发布、每次部署、每次依赖更新时保证：自身及其代为托管、镜像、构建本作品的任何一方，均不是受限实体。

## 六、贡献

6.1 贡献者保证自己不是指定人员或其近亲属，也不代表受限实体行事。

6.2 违反 6.1 的贡献自始无效，可依 4.2 移除。

## 七、救济

7.1 违约的，维护者有权要求禁令救济、损害赔偿、返还因访问获得的收益及维权费用，以适用法律允许的最大限度为限。各方同意未经许可的访问会造成不可弥补的损害。

7.2 任何条款被认定不可执行的，在法律允许的最大范围内仍然有效，其余条款继续有效。
